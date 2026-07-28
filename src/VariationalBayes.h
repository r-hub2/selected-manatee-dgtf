#ifndef _VARIATIONALBAYES_HPP
#define _VARIATIONALBAYES_HPP

#include <chrono>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <algorithm>
#include <RcppArmadillo.h>
#include "Model.h"
#include "LinkFunc.h"
#include "yjtrans.h"
#include "LinearBayes.h"
#include "SequentialMonteCarlo.h"
#include "StaticParams.h"

#ifdef DGTF_TIMING_HVA
#include <chrono>
#define T_NOW() std::chrono::high_resolution_clock::now()
#define T_US(dt) std::chrono::duration_cast<std::chrono::microseconds>(dt).count()
#endif


// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(RcppArmadillo)]]

namespace VB
{
    class VariationalBayes
    {
    public:
        unsigned int nsample = 1000;
        // unsigned int nthin = 2;
        unsigned int niter = 1000;
        // unsigned int ntotal = 3001;
        unsigned int N = 500; // number of SMC particles
        unsigned int N_sample = 50;

        bool use_discount = false;
        double discount_factor = 0.95;

        bool update_static = true;
        unsigned int m = 1; // number of unknown static parameters
        std::vector<std::string> param_selected = {"W"};

        Prior W_prior, seas_prior, rho_prior, par1_prior, par2_prior;
        bool zintercept_infer = false;
        bool zzcoef_infer = false;

        arma::vec W_stored;    // nsample x 1
        arma::mat seas_stored;  // period x nsample
        arma::vec rho_stored;  // nsample x 1
        arma::vec par1_stored; // nsample x 1
        arma::vec par2_stored; // nsample x 1
        arma::mat psi_stored; // (nT + 1) x nsample
        arma::mat z_stored; // (nT + 1) x nsample
        arma::mat prob_stored; // (nT + 1) x nsample
        arma::vec zintercept_stored; // nsample x 1
        arma::vec zzcoef_stored; // nsample x 1

        VariationalBayes(const Model &model, const Rcpp::List &vb_opts)
        {
            Rcpp::List opts = vb_opts;

            nsample = 1000;
            if (opts.containsElementNamed("nsample"))
            {
                nsample = Rcpp::as<unsigned int>(opts["nsample"]);
            }

            // nthin = 1000;
            // if (opts.containsElementNamed("nthin"))
            // {
            //     nthin = Rcpp::as<unsigned int>(opts["nthin"]);
            // }

            niter = 5000;
            if (opts.containsElementNamed("niter"))
            {
                niter = Rcpp::as<unsigned int>(opts["niter"]);
            }
            else if (opts.containsElementNamed("nburnin"))
            {
                niter = Rcpp::as<unsigned int>(opts["nburnin"]);
            }

            // ntotal = niter + nthin * nsample + 1;

            N = 500;
            if (opts.containsElementNamed("num_particle"))
            {
                N = Rcpp::as<unsigned int>(opts["num_particle"]);
            }

            N_sample = 50;
            if (opts.containsElementNamed("num_particle_sample"))
            {
                N_sample = Rcpp::as<unsigned int>(opts["num_particle_sample"]);
            }

            use_discount = false;
            if (opts.containsElementNamed("use_discount"))
            {
                use_discount = Rcpp::as<bool>(opts["use_discount"]);
            }

            discount_factor = 0.95;
            if (opts.containsElementNamed("discount_factor"))
            {
                discount_factor = Rcpp::as<double>(opts["discount_factor"]);
            }

            Static::init_prior(
                param_selected, m,
                W_prior, seas_prior, rho_prior,
                par1_prior, par2_prior,
                zintercept_infer, zzcoef_infer,
                opts, model);

            update_static = false;
            if (m > 0)
            {
                update_static = true;
            }

            if (W_prior.infer)
            {
                W_stored.set_size(nsample);
                W_stored.zeros();
            }

            if (seas_prior.infer)
            {
                seas_stored.set_size(model.seas.period, nsample);
                seas_stored.zeros();
            }

            if (rho_prior.infer)
            {
                rho_stored.set_size(nsample);
                rho_stored.zeros();
            }

            if (par1_prior.infer)
            {
                par1_stored.set_size(nsample);
                par1_stored.zeros();
            }

            if (par2_prior.infer)
            {
                par2_stored.set_size(nsample);
                par2_stored.zeros();
            }

            if (zintercept_infer)
            {
                zintercept_stored.set_size(nsample);
                zintercept_stored.zeros();
            }

            if (zzcoef_infer) {
                zzcoef_stored.set_size(nsample);
                zzcoef_stored.zeros();
            }
            
        }


        static Rcpp::List default_settings()
        {
            Rcpp::List opts = Static::default_settings();
            opts["nsample"] = 1000;
            // opts["nthin"] = 1;
            opts["niter"] = 5000;
            opts["num_particle"] = 500;
            opts["use_discount"] = false;
            opts["discount_factor"] = 0.95;
            return opts;
        }
    };
    /**
     * @brief Gradient ascent.
     *
     */
    class HybridParams
    {
    public:
        HybridParams() : change(par_change)
        {
            learning_rate = 0.01;
            eps_step_size = 1.e-6;
        }

        HybridParams(
            const unsigned int &m_in,
            const double &learn_rate = 0.01,
            const double &eps_size = 1.e-6) : change(par_change)
        {
            init(m_in, learn_rate, eps_size);
        }

        HybridParams(
            const arma::vec &curEg2_init,     // m x 1
            const arma::vec &curEdelta2_init, // m x 1
            const double &learn_rate = 0.01,
            const double &eps_size = 1.e-6) : change(par_change)
        {
            init(curEg2_init, curEdelta2_init, learn_rate, eps_size);
        }

        void init(
            const unsigned int &m_in,
            const double &learn_rate = 0.01,
            const double &eps_size = 1.e-6)
        {
            m = m_in;
            learning_rate = learn_rate;
            eps_step_size = eps_size;

            par_change.set_size(m);
            par_change.zeros();

            curEg2 = par_change;
            curEdelta2 = par_change;
            rho = par_change;
        }

        void init(
            const arma::vec &curEg2_init,     // m x 1
            const arma::vec &curEdelta2_init, // m x 1
            const double &learn_rate = 0.01,
            const double &eps_size = 1.e-6)
        {
            m = curEg2_init.n_elem;
            learning_rate = learn_rate;
            eps_step_size = eps_size;

            curEg2 = curEg2_init;
            curEdelta2 = curEdelta2_init;

            par_change.set_size(m);
            par_change.zeros();
            rho = par_change;
        }

        void update_grad(const arma::vec &dYJinv_dVecPar) // Checked. OK.
        {
            curEg2 *= (1.0 - learning_rate);  // In-place multiply
            curEg2 += learning_rate * arma::square(dYJinv_dVecPar);  // In-place add

            rho = arma::sqrt((curEdelta2 + eps_step_size) / (curEg2 + eps_step_size));
            par_change = rho % dYJinv_dVecPar;

            curEdelta2 *= (1.0 - learning_rate); // In-place multiply
            curEdelta2 += learning_rate * arma::square(par_change); // In-place add

            #ifdef DGTF_DO_BOUND_CHECK
                bound_check<arma::vec>(curEg2, "curEg2");
                bound_check<arma::vec>(par_change, "update_grad: par_change");
                bound_check<arma::vec>(curEdelta2, "update_grad: curEdelta2");
            #endif
            
            return;
        }

        const arma::vec &change;

    private:
        arma::vec curEg2;     // m x 1
        arma::vec curEdelta2; // m x 1
        arma::vec par_change;

        arma::vec rho; // m x 1, step size

        double learning_rate = 0.01;
        double eps_step_size = 1.e-6;
        unsigned int m;
    };

    class Hybrid : public VariationalBayes
    {
    public:
        arma::cube Theta_stored; // nP x (nT + 1) x nsample
        std::string fsys;

        Hybrid(const Model &model_in, const Rcpp::List &hvb_opts) : VariationalBayes(model_in, hvb_opts)
        {
            Rcpp::List opts = hvb_opts;

            learning_rate = 0.01;
            if (opts.containsElementNamed("learning_rate"))
            {
                learning_rate = Rcpp::as<double>(opts["learning_rate"]);
            }

            eps_step_size = 1.e-6;
            if (opts.containsElementNamed("eps_step_size"))
            {
                eps_step_size = Rcpp::as<double>(opts["eps_step_size"]);
            }

            k = 1;
            if (opts.containsElementNamed("k"))
            {
                k = Rcpp::as<unsigned int>(opts["k"]);
            }
            k = std::min(k, m);

            marglike_stored.set_size(niter);
            marglike_stored.zeros();

            // condlike_stored.set_size(niter);
            // condlike_stored.zeros();

            // grad_stored.set_size(niter);
            // grad_stored.zeros();

            gamma.set_size(m);
            gamma.ones();
            grad_tau.init(m, learning_rate, eps_step_size);

            mu.set_size(m);
            mu.zeros();
            grad_mu.init(m, learning_rate, eps_step_size);

            d.set_size(m);
            d.ones();
            grad_d.init(m, learning_rate, eps_step_size);

            eps = d;

            eta = Static::init_eta(param_selected, model_in, update_static); // Checked. OK.
            eta_tilde = Static::eta2tilde(
                eta, param_selected, W_prior.name, 
                par1_prior.name, model_in.dobs.name,
                model_in.seas.period, model_in.seas.in_state);

            nu = tYJ(eta_tilde, gamma);

            // Initialize variational mean to match model's initial parameter values.
            // Without this, mu=0 maps to W=exp(0)=1 for half-t/gamma priors,
            // which is far from typical W values (~0.01). The resulting huge
            // initial gradient can trigger a collapse: W shrinks → particles
            // collapse → res≈0 → gradient≈-T/2 → W shrinks further → underflow.
            // Priors with strong barriers (inverse-gamma) survive this; weak
            // priors (half-t, half-Cauchy) do not.
            mu = nu;

            B.set_size(m, k);
            B.zeros();
            grad_vecB.init(m * k, learning_rate, eps_step_size);

            xi.set_size(k);
            xi.zeros();
        }


        static Rcpp::List default_settings()
        {
            Rcpp::List opts = VariationalBayes::default_settings();
            opts["learning_rate"] = 0.01;
            opts["eps_step_size"] = 1.e-6;
            opts["k"] = 1;
            return opts;
        }

       
        void infer(
            Model &model,
            const arma::vec &y,
            const bool &verbose = VERBOSE)
        {
            std::map<std::string, SysEq::Evolution> sys_list = SysEq::sys_list;
            fsys = model.fsys;
            const unsigned int nT = y.n_elem - 1;

            // Iterative transfer function (sys_nbinom): par2 is the lag
            // order r, an integer that sets the state dimension r+1.
            // Inferring it would change nP mid-fit. par1 = kappa is fine
            // to infer - its gradient is computed below by borrowing the
            // equivalent untruncated sliding form (nL = nT).
            if (sys_list[model.fsys] == SysEq::Evolution::nbinom &&
                par2_prior.infer)
            {
                throw std::invalid_argument(
                    "VB::Hybrid::infer: priors on lag par2 (the lag order r) "
                    "are not supported with sys_nbinom() under HVA. The state "
                    "dimension r+1 is fixed at construction; inferring r "
                    "would change it mid-fit.");
            }

            // Likewise, full-rank system noise adds Gaussian innovations
            // to every state position. Under the iterative convention,
            // positions 1..r are deterministic f-values; injecting noise
            // there is inconsistent with the model.
            if (sys_list[model.fsys] == SysEq::Evolution::nbinom &&
                model.derr.full_rank)
            {
                throw std::invalid_argument(
                    "VB::Hybrid::infer: sys_nbinom() requires univariate "
                    "system noise (full_rank = FALSE); the iterative state "
                    "carries deterministic past-f values that cannot be "
                    "perturbed independently of psi.");
            }

            model.zero.z.set_size(y.n_elem);
            model.zero.z.ones();
            model.zero.prob = model.zero.z;

            arma::mat Theta(model.nP, y.n_elem, arma::fill::zeros); // nP x (nT + 1)
            arma::cube Theta_all(model.nP, N, y.n_elem);
            arma::vec z(y.n_elem, arma::fill::ones);

            auto start = std::chrono::high_resolution_clock::now();
            for (unsigned int b = 0; b < niter; b++)
            {
                Rcpp::checkUserInterrupt();

                // TFS sampler
                // ------------------
                // MUST set initial_resample_all = true (MCS smoothing) and final_resample_by_weights = false (reduce degeneracy) to make this algorithm work.

                Theta.zeros();
                z.ones();
                double marg_loglik = SMC::SequentialMonteCarlo::auxiliary_filter0(
                    Theta, Theta_all, z, model, y, N, 
                    true, false, use_discount, discount_factor);
                // arma::mat Theta = arma::mean(Theta_tmp, 1); // nP x (nT + 1)

                marglike_stored.at(b) = marg_loglik;

                if (model.zero.inflated)
                {
                    model.zero.prob = z; // (nT + 1) x 1
                    // Vectorized Bernoulli using Armadillo RNG (no R::runif loop)
                    arma::vec u = arma::randu<arma::vec>(model.zero.z.n_elem);
                    model.zero.z = arma::conv_to<arma::vec>::from(u < model.zero.prob);
                }
                // ------------------

                // Compute gradient
                arma::mat dFphi_grad;
                arma::vec dloglik_dlag(2, arma::fill::zeros);
                const bool iterative_ft =
                    (sys_list[model.fsys] == SysEq::Evolution::nbinom);
                // par2 (= r) is structural under iterative and refused at
                // the entry guard; only par1 (= kappa) reaches here. For
                // sliding, both par1 and par2 are continuous knobs.
                const bool need_lag_grad =
                    par1_prior.infer || (par2_prior.infer && !iterative_ft);
                if (need_lag_grad)
                {
                    // For the iterative form (sys_nbinom), the kappa-gradient
                    // of f_t is computed by borrowing the untruncated sliding
                    // sum: the iterative recursion and the sliding sum with
                    // nL = nT compute identical f_t, so their gradients
                    // w.r.t. logit(kappa) coincide. We use the sliding
                    // form here because its gradient is parallel-safe over
                    // t (no sequential recursion). LagDist::get_Fphi_grad
                    // returns derivatives w.r.t. logit(kappa) (par1) and
                    // log(r) (par2) - the same unconstrained parameterisa-
                    // tion that downstream Static::dlogJoint_deta expects.
                    const unsigned int grad_nL =
                        iterative_ft ? nT : model.dlag.nL;
                    dFphi_grad = LagDist::get_Fphi_grad(
                        grad_nL, model.dlag.name,
                        model.dlag.par1, model.dlag.par2);
                }

                arma::vec lambda(y.n_elem, arma::fill::zeros);
                double dpar1_acc = 0.0, dpar2_acc = 0.0;
                // Precompute seasonal X^T * val once per iteration (val changes across b, not within)
                arma::vec seas_xt;
                const bool add_seas =
                    (model.seas.period > 0) && (!model.seas.in_state) &&
                    (!model.seas.X.is_empty()) && (!model.seas.val.is_empty());
                if (add_seas) {
                    seas_xt = model.seas.X.t() * model.seas.val; // (nT+1) x 1
                }

                // Fast access to Fphi and (optional) its gradients
                const double *dF1_ptr = need_lag_grad ? dFphi_grad.colptr(0) : nullptr;
                const double *dF2_ptr = need_lag_grad ? dFphi_grad.colptr(1) : nullptr;
                const double *Fphi_ptr = model.dlag.Fphi.memptr();
                const unsigned int iter_grad_nelem_max = iterative_ft ? nT : 0u;
                #ifdef DGTF_USE_OPENMP
                #pragma omp parallel for schedule(static) reduction(+ : dpar1_acc, dpar2_acc)
                #endif
                for (unsigned int t = 1; t < y.n_elem; t++)
                {
                    if (model.zero.z.at(t) < EPS)
                        continue;

                    const double *theta_col = Theta.colptr(t);
                    double ft = 0.0, dpar1_sum = 0.0, dpar2_sum = 0.0;

                    if (iterative_ft)
                    {
                        // Iterative state: theta[t] = (psi[t+1], f[t],
                        // f[t-1], ..., f[t+1-r]). The transfer-function
                        // value f[t] is stored directly at position 1
                        // by SysEq::func_gt; read it out, no recomputation.
                        ft = theta_col[1];

                        if (need_lag_grad)
                        {
                            // Borrow the untruncated sliding gradient
                            // (nL = nT): the iterative f_t equals the
                            // sliding sum with no truncation, so their
                            // partial derivatives w.r.t. logit(kappa)
                            // coincide. Under iterative state convention
                            // psi[t-l] lives at Theta.at(0, t-1-l), not
                            // at theta_col[l] - row 0 of column t-1-l.
                            const unsigned int nelem =
                                std::min(t, iter_grad_nelem_max);
                            for (unsigned int l = 0; l < nelem; ++l)
                            {
                                const double ylag = y.at(t - 1 - l);
                                if (std::abs(ylag) > 0.0)
                                {
                                    const double psi_lag = Theta.at(0, t - 1 - l);
                                    const double hpsi = GainFunc::psi2hpsi(psi_lag, model.fgain);
                                    dpar1_sum += dF1_ptr[l] * hpsi * ylag;
                                }
                            }
                            // dpar2_sum stays 0: par2 (= r) is structural
                            // under iterative and refused at the entry guard.
                        }
                    }
                    else
                    {
                        const unsigned int nelem = std::min(t, model.dlag.nL);

                        // Manual loop unrolling for small nelem
                        unsigned int i = 0;
                        for (; i + 4 <= nelem; i += 4)
                        {
                            // Process 4 at a time
                            for (int j = 0; j < 4; ++j)
                            {
                                const unsigned int idx = i + j;
                                const double ylag = y.at(t - 1 - idx);
                                if (std::abs(ylag) > 0.0)
                                {
                                    const double psi = theta_col[idx];
                                    const double hpsi = GainFunc::psi2hpsi(psi, model.fgain);
                                    const double common = hpsi * ylag;

                                    ft += Fphi_ptr[idx] * common;
                                    if (need_lag_grad)
                                    {
                                        dpar1_sum += dF1_ptr[idx] * common;
                                        dpar2_sum += dF2_ptr[idx] * common;
                                    }
                                }
                            }
                        }

                        // Handle remainder
                        for (; i < nelem; ++i)
                        {
                            const double ylag = y.at(t - 1 - i);
                            if (std::abs(ylag) > 0.0)
                            {
                                const double psi = theta_col[i];
                                const double hpsi = GainFunc::psi2hpsi(psi, model.fgain);
                                const double common = hpsi * ylag;

                                ft += Fphi_ptr[i] * common;
                                if (need_lag_grad)
                                {
                                    dpar1_sum += dF1_ptr[i] * common;
                                    dpar2_sum += dF2_ptr[i] * common;
                                }
                            }
                        }
                    }

                    if (add_seas)
                    {
                        ft += seas_xt.at(t);
                    }

                    lambda.at(t) = LinkFunc::ft2mu(ft, model.flink);

                    const double dll_deta = Model::dloglik_deta(
                        ft, y.at(t), model.dobs.par2, model.dobs.name, model.flink);

                    dpar1_acc += dll_deta * dpar1_sum;
                    dpar2_acc += dll_deta * dpar2_sum;
                }

                if (par1_prior.infer || par2_prior.infer)
                {
                    dloglik_dlag.at(0) = dpar1_acc;
                    dloglik_dlag.at(1) = dpar2_acc;
                }

                // HVB parameter updates
                if (update_static)
                {
                    arma::vec dlogJoint = Static::dlogJoint_deta(
                        y, Theta, lambda, dloglik_dlag, eta, param_selected,
                        W_prior, par1_prior, par2_prior, rho_prior, seas_prior, model); // Checked. OK.

                    arma::mat SigInv = get_sigma_inv(B, d, k);
                    arma::vec dlogq = dlogq_dtheta(SigInv, nu, eta_tilde, gamma, mu);
                    arma::vec ddiff = dlogJoint - dlogq;

                    // mu
                    arma::vec dyji_dnu = dYJinv_dnu_diag(nu, gamma); // m x 1
                    arma::vec L_mu = dyji_dnu % ddiff;
                    grad_mu.update_grad(L_mu);
                    mu += grad_mu.change;

                    // grad_stored.at(b) += arma::accu(arma::abs(grad_mu.change));

                    if (m > 1)
                    {
                        arma::mat L_B(m, k, arma::fill::zeros);
                        if (k > 1)
                        {
                            L_B = L_mu * xi.t(); // m x k, Appendix B.2 equation (ii)

                            // Enforce lower-triangular (or diagonal) constraint
                            for (unsigned int col = 0; col < k; ++col) {
                                for (unsigned int row = 0; row < col; ++row) { // strict upper
                                    L_B(row, col) = 0.0;
                                }
                            }
                            
                        }
                        else
                        {
                            L_B.col(0) = L_mu * xi[0];
                        }


                        // Non‑owning view over memory (no copy)
                        arma::vec L_B_view(const_cast<double*>(L_B.memptr()), m * k, false, true);
                        grad_vecB.update_grad(L_B_view);

                        // Apply update (reshape view of change)
                        arma::mat B_change2(const_cast<double*>(grad_vecB.change.memptr()), m, k, false, true);
                        if (k == 1) {
                            B.col(0) += B_change2.col(0);
                        } else {
                            for (unsigned int col = 0; col < k; ++col) {
                                for (unsigned int row = 0; row < m; ++row) {
                                    if (row < col) {
                                        B(row, col) = 0.0; // strict upper -> zero
                                    } else {
                                        B(row, col) += B_change2(row, col); // lower/diag -> add
                                    }
                                }
                            }
                        }

                        // grad_stored.at(b) += arma::accu(arma::abs(grad_vecB.change));
                    }

                    // d
                    arma::vec L_d = (eps % dyji_dnu) % ddiff; // m x 1, Section B.2 equation (iii), Ref: `dtheta_dBDelta.m`
                    grad_d.update_grad(L_d);
                    d += grad_d.change;
                    // grad_stored.at(b) += arma::accu(arma::abs(grad_d.change));

                    // tau
                    arma::vec tau = gamma2tau(gamma);
                    arma::vec L_tau = dYJinv_dtau_diag(nu, gamma) % ddiff;
                    grad_tau.update_grad(L_tau);
                    tau += grad_tau.change;
                    tau2gamma(tau, gamma);
                    // grad_stored.at(b) += arma::accu(arma::abs(grad_tau.change));

                    rtheta(nu, eta_tilde, xi, eps, gamma, mu, B, d);

                    // Safety floor for variance-type parameters in unconstrained space.
                    // Even with proper mu initialization, noisy SMC gradients can
                    // occasionally push eta_tilde too far negative, triggering the
                    // collapse spiral (small W → res≈0 → gradient≈-T/2 → smaller W).
                    // Floor at exp(-20) ≈ 2e-9, well below any realistic W but far
                    // above double underflow (~exp(-745)).
                    {
                        std::map<std::string, AVAIL::Dist> dist_list_local = AVAIL::dist_list;
                        std::map<std::string, AVAIL::Param> sp_list = AVAIL::static_param_list;
                        const double eta_tilde_floor = -20.0;
                        unsigned int jj = 0;
                        for (unsigned int kk = 0; kk < param_selected.size(); kk++)
                        {
                            switch (sp_list[param_selected[kk]])
                            {
                            case AVAIL::Param::W:
                            {
                                // For half-t/gamma/halfcauchy: eta_tilde = log(W)
                                // For invgamma: eta_tilde = -log(W), floor is not needed
                                //   (the -beta/W barrier prevents collapse)
                                if (dist_list_local[W_prior.name] != AVAIL::Dist::invgamma)
                                {
                                    if (eta_tilde.at(jj) < eta_tilde_floor)
                                        eta_tilde.at(jj) = eta_tilde_floor;
                                }
                                jj += 1;
                                break;
                            }
                            case AVAIL::Param::seas:
                            {
                                if (!model.seas.in_state)
                                    jj += model.seas.period;
                                break;
                            }
                            default:
                            {
                                jj += 1;
                                break;
                            }
                            }
                        }
                    }

                    eta = Static::tilde2eta(
                        eta_tilde, param_selected, 
                        W_prior.name, par1_prior.name, 
                        model.dlag.name, model.dobs.name,
                        model.seas.period, model.seas.in_state);
                    Static::update_params(model, param_selected, eta);
                } // end update_static

                if (verbose)
                {
                    Rprintf("\rOptimization: %u/%u", b + 1, niter);
                }

            } // HVB SGD Loop

            auto end = std::chrono::high_resolution_clock::now();
            elapsed_opt_us = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
            Rprintf("\nOptimization done in %.2f s.\n", elapsed_opt_us / 1e6);


            if (verbose)
            {
                Rprintf("\n");
            }

            // ------------------------------------------------------------------
            // For discounted runs, the per-time-step innovation variance W_t is
            // computed adaptively by LinearBayes during inference rather than
            // estimated as a single scalar. The forecast routine in
            // src/export.cpp::dgtf_forecast() needs a representative scalar W
            // to drive the leading-psi random walk past t = T; without one it
            // falls back to model.derr.par1 (the initial constructor W), which
            // is generally not what inference used. Compute the last in-sample
            // W_t here using the post-optimization model parameters and stash
            // it as a class member for get_output() to expose.
            //
            // Wt is shaped nP x nP x (T + 1); the leading-component variance
            // Wt(0, 0, T) is the scalar that drives the shift-system random
            // walk on psi.
            // ------------------------------------------------------------------
            if (use_discount)
            {
                LBA::LinearBayes lba(use_discount, discount_factor);
                lba.filter(model, y);
                arma::cube Wt = lba.get_Wt(model, y, discount_factor);
                if (Wt.n_slices > 0 && Wt.n_rows > 0 && Wt.n_cols > 0)
                {
                    W_forecast_ = Wt(0, 0, Wt.n_slices - 1);
                    W_forecast_set_ = true;
                }
            }

            if (model.zero.inflated)
            {
                z_stored.set_size(y.n_elem, nsample);
                prob_stored.set_size(y.n_elem, nsample);
            }

            if (sys_list[model.fsys] == SysEq::Evolution::identity)
            {
                psi_stored.set_size(model.nP, nsample);
                Theta_stored = arma::zeros<arma::cube>(model.nP, y.n_elem, nsample);
            }
            else
            {
                psi_stored.set_size(y.n_elem, nsample);
            }

            psi_stored.zeros();

            std::map<std::string, AVAIL::Param> static_param_list = AVAIL::static_param_list;

            // arma::cube Theta_tmp = arma::zeros<arma::cube>(model.nP, nsample, y.n_elem);
            // arma::mat ztmp(nsample, y.n_elem, arma::fill::ones);
            // arma::mat ptmp(nsample, y.n_elem, arma::fill::randu);
            // double log_cond_marg = SMC::SequentialMonteCarlo::auxiliary_filter0(
            //     Theta_tmp, ztmp, ptmp, model, y, nsample,
            //     true, true, false, 1.);
            
            // if (sys_list[model.fsys] != SysEq::Evolution::identity)
            // {
            //     psi_stored = Theta_tmp.row_as_mat(0); // (nT + 1) x nsample
            // }


            arma::mat eta_tilde = rtheta_batch(gamma, mu, B, d, nsample); // m x nsample

            // Apply the same safety floor as in the SGD loop to sampled eta_tilde values
            {
                std::map<std::string, AVAIL::Dist> dist_list_local = AVAIL::dist_list;
                std::map<std::string, AVAIL::Param> sp_list = AVAIL::static_param_list;
                const double eta_tilde_floor = -20.0;
                unsigned int jj = 0;
                for (unsigned int kk = 0; kk < param_selected.size(); kk++)
                {
                    switch (sp_list[param_selected[kk]])
                    {
                    case AVAIL::Param::W:
                    {
                        if (dist_list_local[W_prior.name] != AVAIL::Dist::invgamma)
                        {
                            eta_tilde.row(jj) = arma::clamp(eta_tilde.row(jj), eta_tilde_floor, UPBND);
                        }
                        jj += 1;
                        break;
                    }
                    case AVAIL::Param::seas:
                    {
                        if (!model.seas.in_state)
                            jj += model.seas.period;
                        break;
                    }
                    default:
                    {
                        jj += 1;
                        break;
                    }
                    }
                }
            } // end eta_tilde floor

            arma::cube Theta_all_sample(model.nP, N_sample, y.n_elem);
            arma::mat Theta_sample(model.nP, y.n_elem, arma::fill::zeros);
            arma::vec z_sample(y.n_elem, arma::fill::ones);

            start = std::chrono::high_resolution_clock::now();
            for (unsigned int i = 0; i < nsample; i++)
            {
                eta = Static::tilde2eta(
                    eta_tilde.col(i), param_selected, 
                    W_prior.name, par1_prior.name, 
                    model.dlag.name, model.dobs.name,
                    model.seas.period, model.seas.in_state);
                
                Static::update_params(model, param_selected, eta);

                Theta_sample.zeros();
                z_sample.ones();
                double log_cond_marg = SMC::SequentialMonteCarlo::auxiliary_filter0(
                    Theta_sample, Theta_all_sample, z_sample, model, y, N_sample,
                    true, false, use_discount, discount_factor);

                if (model.zero.inflated)
                {
                    prob_stored.col(i) = z_sample; // (nT + 1) x 1
                    // Vectorized Bernoulli for storage
                    arma::vec u = arma::randu<arma::vec>(model.zero.z.n_elem);
                    z_stored.col(i) = arma::conv_to<arma::vec>::from(u < prob_stored.col(i));
                }


                if (sys_list[model.fsys] == SysEq::Evolution::identity)
                {
                    psi_stored.col(i) = Theta_sample.col(y.n_elem - 1);
                    Theta_stored.slice(i) = Theta_sample;
                }
                else
                {
                    psi_stored.col(i) = arma::vectorise(Theta_sample.row(0));
                }


                unsigned int idx = 0;
                for (unsigned int k = 0; k < param_selected.size(); k++)
                {
                    /*
                    idx: location of parameter in eta
                    k: location of parameter name in param_selected
                    i: the i-th sample in **_stored
                    */
                    double val = eta.at(idx);
                    switch (static_param_list[param_selected[k]])
                    {
                    case AVAIL::Param::W:
                    {
                        W_stored.at(i) = val;
                        idx += 1;
                        break;
                    }
                    case AVAIL::Param::seas:
                    {
                        if (!model.seas.in_state)
                        {
                            arma::vec seass = eta.subvec(idx, idx + model.seas.period - 1);
                            seas_stored.col(i) = seass;
                            idx += model.seas.period;
                        }
                        break;
                    }
                    case AVAIL::Param::rho:
                    {
                        rho_stored.at(i) = val;
                        idx += 1;
                        break;
                    }
                    case AVAIL::Param::lag_par1:
                    {
                        par1_stored.at(i) = val;
                        idx += 1;
                        break;
                    }
                    case AVAIL::Param::lag_par2:
                    {
                        par2_stored.at(i) = val;
                        idx += 1;
                        break;
                    }
                    case AVAIL::Param::zintercept:
                    {
                        zintercept_stored.at(i) = val;
                        idx += 1;
                        break;
                    }
                    case AVAIL::Param::zzcoef:
                    {
                        zzcoef_stored.at(i) = val;
                        idx += 1;
                        break;
                    }
                    default:
                    {
                        throw std::invalid_argument("VariationalBayes::infer: undefined static parameter " + param_selected[k]);
                    }
                    }
                }

                if (verbose)
                {
                    Rprintf("\rSampling: %u/%u", i + 1, nsample);
                }
            } // sampling loop


            end = std::chrono::high_resolution_clock::now();
            elapsed_sample_us = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
            Rprintf("\nSampling done in %.2f s.\n", elapsed_sample_us / 1e6);

            if (verbose)
            {
                Rprintf("\n");
            }
        } // end infer

        Rcpp::List get_output()
        {
            Rcpp::List output;
            output["marglik"] = Rcpp::wrap(marglike_stored);
            // output["condlik"] = Rcpp::wrap(condlike_stored);
            // output["grad_norm"] = Rcpp::wrap(grad_stored);

            output["psi_stored"] = Rcpp::wrap(psi_stored);
            arma::vec qprob = {0.025, 0.5, 0.975};
            arma::mat psi_quantile = arma::quantile(psi_stored, qprob, 1);
            output["psi"] = Rcpp::wrap(psi_quantile);

            std::map<std::string, SysEq::Evolution> sys_list = SysEq::sys_list;
            if (sys_list[fsys] == SysEq::Evolution::identity)
            {
                output["Theta"] = Rcpp::wrap(Theta_stored);
            }

            if (!z_stored.is_empty() && !prob_stored.is_empty())
            {
                output["z"] = Rcpp::wrap(arma::vectorise(arma::mean(z_stored, 1)));
                output["prob"] = Rcpp::wrap(prob_stored);
            }

            if (W_prior.infer)
            {
                output["W"] = Rcpp::wrap(W_stored.t());
            }

            if (seas_prior.infer)
            {
                output["seas"] = Rcpp::wrap(seas_stored);
            }

            if (rho_prior.infer)
            {
                output["rho"] = Rcpp::wrap(rho_stored.t());
            }

            if (par1_prior.infer)
            {
                output["par1"] = Rcpp::wrap(par1_stored.t());
            }

            if (par2_prior.infer)
            {
                output["par2"] = Rcpp::wrap(par2_stored.t());
            }

            if (zintercept_infer)
            {
                output["zintercept"] = Rcpp::wrap(zintercept_stored.t());
            }

            if (zzcoef_infer)
            {
                output["zzcoef"] = Rcpp::wrap(zzcoef_stored.t());
            }

            if (W_forecast_set_)
            {
                output["W_forecast"] = Rcpp::wrap(W_forecast_);
            }

            output["inferred"] = Rcpp::wrap(param_selected);
            output["elapsed_opt_us"]    = (double) elapsed_opt_us;
            output["elapsed_sample_us"] = (double) elapsed_sample_us;
            return output;
        }

    private:
        long long elapsed_opt_us = 0;
        long long elapsed_sample_us = 0;

        double learning_rate = 0.01;
        double eps_step_size = 1.e-6;
        unsigned int k = 1; // rank of unknown static parameters.

        // StaticParam W, mu0, delta, kappa, r;
        HybridParams grad_mu, grad_vecB, grad_d, grad_tau;
        arma::vec mu, d, gamma, nu, eps, eta, eta_tilde; // m x 1
        arma::vec marglike_stored; //, condlike_stored, grad_stored; // niter x 1

        double W_forecast_     = 0.0;
        bool W_forecast_set_ = false;

        arma::vec xi;                                    // k x 1
        arma::mat B;                                     // m x k
    }; // class Hybrid
}

#endif
