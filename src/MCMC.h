#ifndef _MCMC_HPP
#define _MCMC_HPP

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <algorithm>
#include <RcppArmadillo.h>
#include <pg.h>
#include "Model.h"
#include "LinkFunc.h"
#include "LinearBayes.h"
#include "StaticParams.h"

// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(RcppArmadillo,pg)]]


// ============================================================================
// HMC Adaptation Structs (Armadillo-based)
// ============================================================================

/**
 * @brief Diagonal mass matrix adapter using Welford's online variance algorithm
 *
 * Tracks running mean and variance of unconstrained parameters during burnin.
 * The estimated mass matrix diagonal approximates the marginal precision of each 
 * parameter, improving HMC sampling efficiency.
 */
struct MassAdapter_arma {
    arma::vec mean;
    arma::vec M2;           // Sum of squared deviations (Welford)
    unsigned int count = 0;

    MassAdapter_arma() = default;

    void init(const arma::vec &initial_params) {
        mean = initial_params;
        M2 = arma::zeros(initial_params.n_elem);
        count = 0;
    }

    void update(const arma::vec &log_params) {
        count++;
        arma::vec delta = log_params - mean;
        mean += delta / static_cast<double>(count);
        M2 += delta % (log_params - mean);
    }

    arma::vec get_mass_diag() const {
        if (count < 10) return arma::ones(mean.n_elem);
        arma::vec var = M2 / static_cast<double>(count - 1);
        // Mass = 1/variance, clamped to [0.1, 100]
        return arma::clamp(1.0 / var, 0.1, 100.0);
    }
};


/**
 * @brief Nesterov dual averaging for HMC step size adaptation
 *
 * Implements the dual averaging scheme from Hoffman & Gelman (2014, JMLR)
 * to automatically tune the leapfrog step size toward a target acceptance rate.
 * 
 * @note Hoffman & Gelman (2014): https://jmlr.org/papers/v15/hoffman14a.html
 */
struct DualAveraging_arma {
    double target_accept = 0.75;
    double mu_da;           // log(10 * eps0)
    double log_eps_bar;     // Smoothed log step size
    double log_eps;         // Current log step size
    double h_bar = 0.0;
    double gamma_da = 0.1;
    double t0_da = 10.0;
    double kappa_da = 0.75;
    unsigned int adapt_count = 0;
    unsigned int min_leaps = 3;
    unsigned int max_leaps = 128;

    DualAveraging_arma() : mu_da(std::log(10 * 0.01)),
                           log_eps_bar(std::log(0.01)),
                           log_eps(std::log(0.01)) {}

    DualAveraging_arma(double eps_init, double target = 0.75)
        : target_accept(target),
          mu_da(std::log(10 * eps_init)),
          log_eps_bar(std::log(eps_init)),
          log_eps(std::log(eps_init)),
          h_bar(0.0), gamma_da(0.1), t0_da(10.0), kappa_da(0.75),
          adapt_count(0), min_leaps(3), max_leaps(128) {}

    double update_step_size(double accept_prob) {
        adapt_count++;
        double t = static_cast<double>(adapt_count);
        h_bar = (1.0 - 1.0 / (t + t0_da)) * h_bar +
                (1.0 / (t + t0_da)) * (target_accept - accept_prob);
        log_eps = mu_da - (std::sqrt(t) / gamma_da) * h_bar;

        double eps = std::min(std::max(std::exp(log_eps), 1e-3), 1.0);
        double w = std::pow(t, -kappa_da);
        log_eps_bar = w * std::log(eps) + (1.0 - w) * log_eps_bar;
        return eps;
    }

    void finalize(double &step_size, unsigned int &nleapfrog, double T_target) const {
        step_size = std::min(std::max(std::exp(log_eps_bar), 1e-3), 1.0);
        unsigned int nlf = static_cast<unsigned int>(std::lround(T_target / step_size));
        nleapfrog = std::max(min_leaps, std::min(max_leaps, nlf));
    }
};


/**
 * @brief Diagnostics storage for 1-D HMC chains
 */
struct HMCDiagnostics_arma {
    double accept_count = 0.0;
    arma::vec energy_diff;
    arma::vec grad_norm;
    arma::vec step_size_stored;
    arma::vec nleapfrog_stored;

    HMCDiagnostics_arma() = default;

    HMCDiagnostics_arma(unsigned int ntotal, unsigned int nburnin) {
        energy_diff = arma::zeros(ntotal);
        grad_norm = arma::zeros(ntotal);
        step_size_stored = arma::zeros(nburnin + 1);
        nleapfrog_stored = arma::zeros(nburnin + 1);
    }

    Rcpp::List to_list() const {
        Rcpp::List out;
        out["energy_diff"] = Rcpp::wrap(energy_diff);
        out["grad_norm"] = Rcpp::wrap(grad_norm);
        if (step_size_stored.n_elem > 0 && arma::accu(arma::abs(step_size_stored)) > EPS)
        {
            out["leapfrog_step_size"] = Rcpp::wrap(step_size_stored);
            out["nleapfrog"] = Rcpp::wrap(nleapfrog_stored);
        }
        return out;
    }
};


namespace MCMC
{
    class Leapfrog
    {
    public:
        /**
         * @brief Calculate gradient of the potential energery. The potential energy is the negative log joint probability of the model.
         * 
         * @note Duane et al. (1987): https://www.sciencedirect.com/science/article/abs/pii/037026938791197X
         * @note Neal (2011): https://arxiv.org/pdf/1206.1901
         * @note Betancourt (2018): https://arxiv.org/abs/1701.02434
         * @note How people introduce their use of HMC - Patel et al. (2021): https://arxiv.org/pdf/2110.08363
         * @note Practical guide to HMC: https://bjlkeng.io/posts/hamiltonian-monte-carlo/
         * 
         * @param q 
         * @param model 
         * @param param_selected 
         * @param W_prior 
         * @param seas_prior 
         * @param rho_prior 
         * @param par1_prior 
         * @param par2_prior 
         * @param zintercept_infer 
         * @param zzcoef_infer 
         * @return arma::vec 
         */
        static arma::vec grad_U(
            const Model &model,
            const arma::vec &params, 
            const arma::vec &y,
            const arma::vec &hpsi,
            const arma::mat &Theta,
            const std::vector<std::string> &param_selected,
            const Prior &W_prior, 
            const Prior &seas_prior, 
            const Prior &rho_prior, 
            const Prior &par1_prior, 
            const Prior &par2_prior, 
            const bool &zintercept_infer, 
            const bool &zzcoef_infer,
            const arma::vec &psi_nc = arma::vec(),
            const bool w_noncentered = false
        )
        {
            arma::vec lambda(y.n_elem, arma::fill::zeros);
            for (unsigned int t = 1; t < y.n_elem; t++)
            {
                double eta = TransFunc::transfer_sliding(t, model.dlag.nL, y, model.dlag.Fphi, hpsi);
                if (model.seas.period > 0)
                {
                    eta += arma::dot(model.seas.X.col(t), model.seas.val);
                }
                lambda.at(t) = LinkFunc::ft2mu(eta, model.flink);
            }

            arma::vec dloglik_dlag = Model::dloglik_dlag(
                y, hpsi, model.dlag.nL, model.dlag.name, model.dlag.par1, model.dlag.par2,
                model.dobs, model.seas, model.zero, model.flink);

            // Compute non-centered W likelihood gradient if requested
            double dloglik_dW_nc_val = std::numeric_limits<double>::quiet_NaN();
            if (w_noncentered && W_prior.infer)
            {
                // v_nc[s] = h'(psi_s) * (-psi_s / 2)
                // deta_t/dWtilde = transfer_sliding(t, nL, y, Fphi, v_nc)
                arma::vec dhpsi = GainFunc::psi2dhpsi<arma::vec>(psi_nc, model.fgain);
                arma::vec v_nc = dhpsi % (-psi_nc * 0.5);

                dloglik_dW_nc_val = 0.0;
                double nc_cnt = 0.0;
                for (unsigned int t = 1; t < y.n_elem; t++)
                {
                    if (model.zero.inflated && model.zero.z.at(t) < EPS) continue;
                    nc_cnt += 1.0;

                    double eta_t = TransFunc::transfer_sliding(t, model.dlag.nL, y, model.dlag.Fphi, hpsi);
                    if (model.seas.period > 0)
                    {
                        eta_t += arma::dot(model.seas.X.col(t), model.seas.val);
                    }
                    double dll_deta = Model::dloglik_deta(
                        eta_t, y.at(t), model.dobs.par2, model.dobs.name, model.flink);

                    double deta_dWtilde = TransFunc::transfer_sliding(
                        t, model.dlag.nL, y, model.dlag.Fphi, v_nc);

                    dloglik_dW_nc_val += dll_deta * deta_dWtilde;
                }
                dloglik_dW_nc_val += 0.5 * nc_cnt;
            }

            arma::vec gd_U = Static::dlogJoint_deta(
                y, Theta, lambda, dloglik_dlag, params, param_selected,
                W_prior, par1_prior, par2_prior, rho_prior, seas_prior, model,
                dloglik_dW_nc_val);
            
            // negate it because U = -log(p) and what we calculate above is the gradient of log(p)
            gd_U *= -1.;

            return gd_U;
        } // grad_U
    }; // class Leapfrog

    class Posterior
    {
    public:
        /**
         * @brief Gibbs sampler for zero-inflated indicator z[t]
         * 
         * @param model 
         * @param y 
         * @param wt 
         */
        static void update_zt(
            Model &model, 
            const arma::vec &y, // (nT + 1) x 1
            const arma::vec &wt
        ) // (nT + 1) x 1
        {
            // lambda: (nT + 1) x 1, conditional expectation of y[t] if z[t] = 1
            arma::vec lambda = model.wt2lambda(y, wt, model.seas.period, model.seas.X, model.seas.val);

            arma::vec p01(y.n_elem, arma::fill::zeros); // p(z[t] = 1 | z[t-1] = 0, gamma)
            p01.at(0) = logistic(model.zero.intercept);
            arma::vec p11(y.n_elem, arma::fill::zeros); // p(z[t] = 1 | z[t-1] = 1, gamma)
            p11.at(0) = logistic(model.zero.intercept + model.zero.coef);

            arma::vec prob_filter(y.n_elem, arma::fill::zeros); // p(z[t] = 1 | y[1:t], gamma)
            prob_filter.at(0) = p01.at(0);
            for (unsigned int t = 1; t < y.n_elem; t++)
            {
                double p0 = model.zero.intercept;
                double p1 = model.zero.intercept + model.zero.coef;
                if (!model.zero.X.is_empty())
                {
                    double val = arma::dot(model.zero.X.col(t), model.zero.beta);
                    p0 += val;
                    p1 += val;
                }
                p01.at(t) = logistic(p0); // p(z[t] = 1 | z[t-1] = 0, gamma)
                p11.at(t) = logistic(p1); // p(z[t] = 1 | z[t-1] = 1, gamma)

                if (y.at(t) > EPS)
                {
                    // If y[t] > 0
                    // We must have p(z[t] = 1 | y[t] > 0) = 1
                    prob_filter.at(t) = 1.;
                }
                else
                {
                    double prob_yzero = ObsDist::loglike(
                        0., model.dobs.name, lambda.at(t), model.dobs.par2, false);

                    double pp1 = prob_filter.at(t - 1) * p11.at(t); // p(z[t-1] = 1) * p(z[t] = 1 | z[t-1] = 1)
                    double pp0 = prob_filter.at(t - 1) * std::abs(1. - p11.at(t)); // p(z[t-1] = 1) * p(z[t] = 0 | z[t-1] = 1)

                    pp1 += std::abs(1. - prob_filter.at(t - 1)) * p01.at(t); // p(z[t-1] = 0) * p(z[t] = 1 | z[t-1] = 0)
                    pp0 += std::abs(1. - prob_filter.at(t - 1)) * std::abs(1. - p01.at(t)); // p(z[t-1] = 0) * p(z[t] = 0 | z[t-1] = 0)

                    pp1 *= prob_yzero; // p(y[t] = 0 | z[t] = 1)
                    prob_filter.at(t) = pp1 / (pp1 + pp0 + EPS);
                }
            } // Forward filtering loop

            model.zero.z.at(y.n_elem - 1) = (R::runif(0., 1.) < prob_filter.at(y.n_elem - 1)) ? 1. : 0.;
            for (unsigned int t = y.n_elem - 2; t > 0; t--)
            {
                if (y.at(t) > EPS)
                {
                    model.zero.z.at(t) = 1.;
                }
                else
                {
                    double p1 = model.zero.z.at(t + 1) > EPS ? p11.at(t + 1) : (1. - p11.at(t + 1)); // p(z[t+1] | z[t] = 1)
                    double prob_backward1 = prob_filter.at(t) * std::abs(p1); // p(z[t] = 1 | y[t] = 0) * p(z[t+1] | z[t] = 1)
                    double p0 = model.zero.z.at(t + 1) > EPS ? p01.at(t + 1) : (1. - p01.at(t + 1)); // p(z[t+1] | z[t] = 0)
                    double prob_backward0 = std::abs(1. - prob_filter.at(t)) * std::abs(p0); // p(z[t] = 0 | y[t] = 0) * p(z[t+1] | z[t] = 0)

                    prob_backward1 = prob_backward1 / (prob_backward1 + prob_backward0 + EPS);
                    model.zero.z.at(t) = (R::runif(0., 1.) < prob_backward1) ? 1. : 0.;
                }
            }

        } // update_zt

        static void update_wt( // Checked. OK.
            arma::vec &wt,     // (nT + 1) x 1
            arma::vec &wt_accept,
            ApproxDisturbance &approx_dlm,
            const arma::vec &y, // (nT + 1) x 1
            Model &model,
            const double &mh_sd = 1.0
        )
        {
            arma::vec ft(y.n_elem, arma::fill::zeros);
            double prior_sd = std::sqrt(model.derr.par1);

            arma::vec lam_old = model.wt2lambda(y, wt, model.seas.period, model.seas.X, model.seas.val);

            for (unsigned int t = 1; t < y.n_elem; t++)
            {
                double wt_old = wt.at(t);

                double logp_old = 0.;
                for (unsigned int i = t; i < y.n_elem; i++)
                {
                    if (!(model.zero.inflated && (model.zero.z.at(i) < EPS)))
                    {
                        // For zero-inflated model, only the y[t] that is not missing taken into account
                        logp_old += ObsDist::loglike(y.at(i), model.dobs.name, lam_old.at(i), model.dobs.par2, true);
                    }
                } // Checked. OK.

                logp_old += R::dnorm4(wt_old, 0., prior_sd, true);

                /*
                Metropolis-Hastings
                */
                approx_dlm.update_by_wt(y, wt);
                arma::vec eta_hat = approx_dlm.get_eta_approx(model.seas); // nT x 1, f0, Fn and psi is updated
                arma::vec lambda_hat = LinkFunc::ft2mu<arma::vec>(eta_hat, model.flink); // nT x 1
                arma::vec Vt_hat = ApproxDisturbance::func_Vt_approx(
                    lambda_hat, model.dobs, model.flink); // nT x 1

                arma::mat Fn = approx_dlm.get_Fn(); // nT x nT
                arma::vec Fnt = Fn.col(t - 1); // nT x 1
                arma::vec Fnt2 = Fnt % Fnt;

                arma::vec tmp = Fnt2 / Vt_hat; // nT x 1, element-wise division
                if (model.zero.inflated)
                {
                    tmp %= model.zero.z.subvec(1, tmp.n_elem); // If y[t] is missing (z[t] = 0), F[t]*F[t]'/V[t] is removed from the posterior variance of the proposal
                }
                double mh_prec = arma::accu(tmp) + EPS8;
                // mh_prec = std::abs(mh_prec) + 1. / w0_prior.par2 + EPS;

                double Bs = 1. / mh_prec;
                double Btmp = std::sqrt(Bs);
                // double Btmp = prior_sd;
                Btmp *= mh_sd;
                // Btmp = std::min(Btmp, 10.);

                double wt_new = R::rnorm(wt_old, Btmp); // Sample from MH proposal
                // bound_check(wt_new, "Posterior::update_wt: wt_new");
                /*
                Metropolis-Hastings
                */

                wt.at(t) = wt_new;
                arma::vec lam = model.wt2lambda(y, wt, model.seas.period, model.seas.X, model.seas.val); // Checked. OK.

                double logp_new = 0.;
                for (unsigned int i = t; i < y.n_elem; i++)
                {
                    if (!(model.zero.inflated && (model.zero.z.at(i) == 0)))
                    {
                        logp_new += ObsDist::loglike(y.at(i), model.dobs.name, lam.at(i), model.dobs.par2, true);
                    }
                } // Checked. OK.

                logp_new += R::dnorm4(wt_new, 0., prior_sd, true); // prior

                double logratio = logp_new - logp_old;
                // logratio += logq_old - logq_new;
                logratio = std::min(0., logratio);
                double logps = 0.;
                if (std::log(R::runif(0., 1.)) < logratio)
                {
                    // accept
                    logps = logp_new;
                    wt_accept.at(t) += 1.;
                    lam_old = lam;
                }
                else
                {
                    // reject
                    logps = logp_old;
                    wt.at(t) = wt_old;
                }
            }
        } // func update_wt

        static double update_static_hmc(
            double &energy_diff_out,
            double &grad_norm_out,
            Model &model,
            const arma::vec &y,
            arma::vec &psi,           // non-const: updated on accept
            arma::vec &wt,            // non-const: updated on accept
            const std::vector<std::string> &param_selected,
            const Prior &W_prior,
            const Prior &seas_prior,
            const Prior &rho_prior,
            const Prior &par1_prior,
            const Prior &par2_prior,
            const double &step_size,
            const arma::vec &mass_diag,
            const bool zintercept_infer = false,
            const bool zzcoef_infer = false,
            const unsigned int &L = 10
        )
        {
            std::map<std::string, AVAIL::Dist> obs_list = ObsDist::obs_list;
            const bool w_nc = W_prior.infer; // use non-centered W when W is inferred

            // ============================================================
            // Non-centered setup: compute standardized cumulative sums
            //   xi[s] = psi[s] / sqrt(W_old),  s = 0, ..., nT
            // so that psi[s] = sqrt(W) * xi[s] for any W during leapfrog.
            // ============================================================
            double W_entry = model.derr.par1;
            double sqrt_W_entry = std::sqrt(W_entry + EPS);
            arma::vec xi;  // cumulative sum of epsilon, fixed during HMC
            if (w_nc)
            {
                xi = psi / sqrt_W_entry;  // xi[0] = 0
            }

            // Current psi, Theta, hpsi (will be recomputed if W moves)
            arma::vec psi_cur = psi;
            arma::mat Theta(model.nP, y.n_elem, arma::fill::zeros);
            Theta.row(0) = psi_cur.t();
            arma::vec hpsi = GainFunc::psi2hpsi<arma::vec>(psi_cur, model.fgain);

            // Save current state for potential revert
            Model mod_current = model;
            arma::vec psi_save = psi;
            arma::vec wt_save = wt;

            // Mass matrix quantities
            arma::vec inv_mass = 1.0 / mass_diag;
            arma::vec sqrt_mass = arma::sqrt(mass_diag);

            // Calculate log of joint probability at current position
            arma::vec lambda(y.n_elem, arma::fill::zeros);
            for (unsigned int t = 1; t < y.n_elem; t++)
            {
                double eta = TransFunc::transfer_sliding(t, model.dlag.nL, y, model.dlag.Fphi, hpsi);
                if (model.seas.period > 0)
                {
                    eta += arma::dot(model.seas.X.col(t), model.seas.val);
                }
                lambda.at(t) = LinkFunc::ft2mu(eta, model.flink);
            }

            double logp_old = Static::logJoint(
                y, Theta, lambda, W_prior, par1_prior, par2_prior, rho_prior, seas_prior,
                model, w_nc);

            // Potential energy: U = -log p
            double current_U = -logp_old;

            // Map parameters to unconstrained space
            arma::vec params = Static::init_eta(param_selected, model, true);
            arma::vec q = Static::eta2tilde(
                params, param_selected, W_prior.name, par1_prior.name,
                model.dobs.name, model.seas.period, model.seas.in_state);

            // Sample momentum from N(0, M) where M = diag(mass_diag)
            arma::vec p(q.n_elem);
            for (unsigned int i = 0; i < q.n_elem; i++)
            {
                p.at(i) = sqrt_mass.at(i) * R::rnorm(0., 1.);
            }

            // Kinetic energy: K = 0.5 * p' * M^{-1} * p
            double current_K = 0.5 * arma::accu(arma::square(p) % inv_mass);

            // Compute initial gradient
            arma::vec grad_U = Leapfrog::grad_U(
                model, params, y, hpsi, Theta, param_selected,
                W_prior, seas_prior, rho_prior, par1_prior, par2_prior,
                zintercept_infer, zzcoef_infer,
                psi_cur, w_nc);

            grad_norm_out = arma::norm(grad_U);

            // Half step for momentum
            p -= 0.5 * step_size * grad_U;

            for (unsigned int i = 1; i <= L; i++)
            {
                // Full step for position: q += eps * M^{-1} * p
                q += step_size * (inv_mass % p);

                params = Static::tilde2eta(
                    q, param_selected, W_prior.name, par1_prior.name,
                    model.dlag.name, model.dobs.name,
                    model.seas.period, model.seas.in_state);
                Static::update_params(model, param_selected, params);

                // Non-centered: recompute psi, Theta, hpsi from new W
                if (w_nc)
                {
                    double sqrt_W_cur = std::sqrt(model.derr.par1 + EPS);
                    psi_cur = sqrt_W_cur * xi;
                    Theta.row(0) = psi_cur.t();
                    hpsi = GainFunc::psi2hpsi<arma::vec>(psi_cur, model.fgain);
                }

                // Recompute gradient
                grad_U = Leapfrog::grad_U(
                    model, params, y, hpsi, Theta, param_selected,
                    W_prior, seas_prior, rho_prior, par1_prior, par2_prior,
                    zintercept_infer, zzcoef_infer,
                    psi_cur, w_nc);

                // Full step for momentum (except at last step)
                if (i != L)
                {
                    p -= step_size * grad_U;
                }
            }

            // Final half step for momentum
            p -= 0.5 * step_size * grad_U;

            // Negate momentum for reversibility
            p *= -1.;

            // Calculate proposed energy (hpsi and psi_cur already up-to-date)
            for (unsigned int t = 1; t < y.n_elem; t++)
            {
                double eta = TransFunc::transfer_sliding(t, model.dlag.nL, y, model.dlag.Fphi, hpsi);
                if (model.seas.period > 0)
                {
                    eta += arma::dot(model.seas.X.col(t), model.seas.val);
                }
                lambda.at(t) = LinkFunc::ft2mu(eta, model.flink);
            }

            double logp_new = Static::logJoint(
                y, Theta, lambda, W_prior, par1_prior, par2_prior, rho_prior, seas_prior,
                model, w_nc);

            double proposed_U = -logp_new;
            double proposed_K = 0.5 * arma::accu(arma::square(p) % inv_mass);

            // Hamiltonian difference
            double H_proposed = proposed_U + proposed_K;
            double H_current = current_U + current_K;
            energy_diff_out = H_proposed - H_current;

            // Metropolis acceptance
            double accept_prob = 0.0;
            bool accept = false;
            if (!std::isfinite(energy_diff_out) || energy_diff_out > 100.0)
            {
                accept_prob = 0.0;
            }
            else if (energy_diff_out <= 0.0)
            {
                accept_prob = 1.0;
                accept = true;
            }
            else
            {
                accept_prob = std::exp(-energy_diff_out);
                if (R::runif(0., 1.) < accept_prob)
                {
                    accept = true;
                }
            }

            if (accept)
            {
                // Non-centered: update psi and wt to reflect the new W
                if (w_nc)
                {
                    psi = psi_cur;
                    wt.at(0) = 0.0;
                    wt.subvec(1, psi.n_elem - 1) = arma::diff(psi);
                }
            }
            else
            {
                // Revert model to current state
                model = mod_current;
                psi = psi_save;
                wt = wt_save;
            }

            return accept_prob;
        } // func update_static_hmc


        static arma::vec update_pg_psi(
            const Model &model, 
            const arma::vec &y, 
            const arma::vec &psi_old // nP x (nT + 1)
        )
        {
            const unsigned int nP = model.nP;
            const unsigned int nT = y.n_elem - 1;
            const double npop = model.dobs.par2;

            Model normal_dlm = model;
            normal_dlm.dobs.name = "gaussian";
            normal_dlm.flink = "identity";
            normal_dlm.fgain = "identity";

            arma::vec omega(nT + 1, arma::fill::ones);
            arma::mat Theta_old = TransFunc::psi2theta(psi_old, y, model.ftrans, model.fgain, model.dlag);
            for (unsigned int t = 0; t <= nT; t++)
            {
                if (std::abs(model.zero.z.at(t) - 1.) < EPS)
                {
                    // If y[t] is not missing
                    double eta = TransFunc::func_ft(model.ftrans, model.fgain, model.dlag, model.seas, t, Theta_old.col(t), y);
                    omega.at(t) = pg::rpg_scalar_hybrid(y.at(t) + npop, eta);
                }
            }

            arma::vec ktmp = 0.5 * (y - npop) / omega;
            // arma::vec ktmp = y;

            // Forward filtering
            arma::mat mt(nP, nT + 1, arma::fill::zeros);
            arma::cube Ct(nP, nP, nT + 1);
            Ct.slice(0) = 5. * arma::eye<arma::mat>(nP, nP);
            arma::mat at = mt;
            arma::cube Rt = Ct;

            arma::vec Ft = TransFunc::init_Ft(nP, model.ftrans, model.seas.period, model.seas.in_state);
            arma::mat Gt = SysEq::init_Gt(nP, model.dlag, model.fsys, model.seas.period, model.seas.in_state);
            for (unsigned int t = 1; t < nT + 1; t++)
            {
                at.col(t) = SysEq::func_gt(
                    model.fsys, model.fgain, model.dlag,
                    mt.col(t - 1), y.at(t - 1),
                    model.seas.period, model.seas.in_state);
                Rt.slice(t) = LBA::func_Rt(Gt, Ct.slice(t - 1), normal_dlm.derr.par1);

                if (std::abs(model.zero.z.at(t) - 1.) < EPS)
                {
                    // When y[t] is not missing
                    normal_dlm.dobs.par2 = 1. / omega.at(t);
                    double ft_prior = 0.;
                    double qt_prior = 0.;
                    LBA::func_prior_ft(ft_prior, qt_prior, Ft, t, normal_dlm, y, at.col(t), Rt.slice(t));
                    qt_prior += normal_dlm.dobs.par2;

                    double ft_posterior = ktmp.at(t);
                    double qt_posterior = 0.;

                    arma::mat At = LBA::func_At(Rt.slice(t), Ft, qt_prior);
                    mt.col(t) = LBA::func_mt(at.col(t), At, ft_prior, ft_posterior);
                    Ct.slice(t) = LBA::func_Ct(Rt.slice(t), At, qt_prior, qt_posterior);
                }
                else
                {
                    // When y[t] is missing
                    mt.col(t) = at.col(t);
                    Ct.slice(t) = Rt.slice(t);
                }
            }

            arma::mat Theta(nP, nT + 1, arma::fill::zeros);
            // Backward sampling
            {
                arma::mat Ct_chol = arma::chol(arma::symmatu(Ct.slice(nT)));
                Theta.col(nT) = mt.col(nT) + Ct_chol.t() * arma::randn(nP);
            }

            for (unsigned int t = nT - 1; t > 0; t--)
            {
                arma::mat Rt_inv = inverse(Rt.slice(t + 1));
                arma::mat Bt = Ct.slice(t) * Gt.t() * Rt_inv;

                arma::vec ht = mt.col(t) + Bt * (Theta.col(t + 1) - at.col(t + 1));
                arma::mat Ht = Ct.slice(t) - Bt * Gt * Ct.slice(t);
                Ht.diag() += EPS8;

                arma::mat Ht_chol = arma::chol(arma::symmatu(Ht));
                Theta.col(t) = ht + Ht_chol.t() * arma::randn(nP);
            }

            arma::vec psi_new = arma::vectorise(Theta.row(0));
            return psi_new;
        } // update_pg_psi
    }; // class Posterior

    class Disturbance
    {
    public:
        Disturbance(const Model &model, const Rcpp::List &mcmc_settings)
        {
            Rcpp::List opts = mcmc_settings;

            L = 50;
            if (opts.containsElementNamed("leap"))
            {
                L = Rcpp::as<unsigned int>(opts["leap"]);
            }
            // Backward compatibility: also accept the earlier "Leap"
            // and the legacy "L" names.
            else if (opts.containsElementNamed("Leap"))
            {
                L = Rcpp::as<unsigned int>(opts["Leap"]);
            }
            else if (opts.containsElementNamed("L"))
            {
                L = Rcpp::as<unsigned int>(opts["L"]);
            }

            max_lag = 50;
            if (opts.containsElementNamed("max_lag"))
            {
                max_lag = Rcpp::as<unsigned int>(opts["max_lag"]);
            }

            mh_sd = 1.0;
            if (opts.containsElementNamed("mh_sd"))
            {
                mh_sd = Rcpp::as<double>(opts["mh_sd"]);
            }

            nburnin = 100;
            if (opts.containsElementNamed("nburnin"))
            {
                nburnin = Rcpp::as<unsigned int>(opts["nburnin"]);
            }
            nthin = 1;
            if (opts.containsElementNamed("nthin"))
            {
                nthin = Rcpp::as<unsigned int>(opts["nthin"]);
            }

            nsample = 100;
            if (opts.containsElementNamed("nsample"))
            {
                nsample = Rcpp::as<unsigned int>(opts["nsample"]);
            }

            ntotal = nburnin + nthin * nsample + 1;

            // HMC adaptation options
            hmc_dual_averaging = true;
            if (opts.containsElementNamed("hmc_dual_averaging"))
            {
                hmc_dual_averaging = Rcpp::as<bool>(opts["hmc_dual_averaging"]);
            }

            hmc_T_target = 2.0;
            if (opts.containsElementNamed("hmc_T_target"))
            {
                hmc_T_target = Rcpp::as<double>(opts["hmc_T_target"]);
            }

            hmc_step_size_init = 0.005;
            if (opts.containsElementNamed("hmc_step_size_init"))
            {
                hmc_step_size_init = Rcpp::as<double>(opts["hmc_step_size_init"]);
            }
            // Backward compatibility: also accept "epsilon" as initial step size
            if (opts.containsElementNamed("epsilon"))
            {
                arma::vec eps = Rcpp::as<arma::vec>(opts["epsilon"]);
                hmc_step_size_init = eps.at(0);
            }

            hmc_step_size = hmc_step_size_init;

            hmc_target_accept = 0.75;
            if (opts.containsElementNamed("hmc_target_accept"))
            {
                hmc_target_accept = Rcpp::as<double>(opts["hmc_target_accept"]);
            }

            Static::init_prior(
                param_selected, nparam,
                W_prior, seas_prior, rho_prior,
                par1_prior, par2_prior,
                zintercept_infer, zzcoef_infer,
                opts, model);

            if (nparam > 0)
            {
                update_static = true;
                hmc_accept = 0.;

                // Initialize mass matrix to identity
                mass_diag.set_size(nparam);
                mass_diag.ones();
            }
            else
            {
                update_static = false;
            }

            W_stored.set_size(nsample);
            // W_accept = 0.;

            seas_stored.set_size(model.seas.period, nsample);
            // seas_accept = 0.;

            rho_stored.set_size(nsample);
            // rho_accept = 0.;

            par1_stored.set_size(nsample);
            par2_stored.set_size(nsample);
            // lag_accept = 0.;

            zintercept_stored.set_size(nsample);
            zzcoef_stored.set_size(nsample);

            return;
        }

        static Rcpp::List default_settings()
        {
            Rcpp::List opts = Static::default_settings();

            opts["epsilon"] = 0.01;
            opts["leap"] = 10;
            opts["max_lag"] = 50;

            opts["mh_sd"] = 1.0;
            opts["nburnin"] = 100;
            opts["nthin"] = 1;
            opts["nsample"] = 100;

            // HMC adaptation
            opts["hmc_dual_averaging"] = true;
            opts["hmc_step_size_init"] = 0.01;
            opts["hmc_T_target"] = 2.0;
            opts["hmc_target_accept"] = 0.75;

            return opts;
        }

        Rcpp::List get_output()
        {
            arma::vec qprob = {0.025, 0.5, 0.975};
            Rcpp::List output;

            arma::mat psi_quantile = arma::quantile(psi_stored, qprob, 1); // (nT + 1) x 3
            output["psi_stored"] = Rcpp::wrap(psi_stored);
            output["psi"] = Rcpp::wrap(psi_quantile);
            output["wt_accept"] = Rcpp::wrap(wt_accept / ntotal);

            if (!z_stored.is_empty())
            {
                output["z_stored"] = Rcpp::wrap(arma::vectorise(arma::mean(z_stored, 1)));
            }
            
            // arma::mat psi_quantile = arma::quantile(wt_stored, qprob, 1); // (nT + 1) x 3
            // output["psi"] = Rcpp::wrap(psi_quantile);
            // output["wt_accept"] = Rcpp::wrap(W_accept / ntotal);
            // output["log_marg"] = Rcpp::wrap(log_marg_stored.t());

            if (update_static)
            {
                output["infer_W"] = W_prior.infer;
                output["W"] = Rcpp::wrap(W_stored);
                // output["W_accept"] = W_accept / ntotal;

                output["infer_seas"] = seas_prior.infer;
                output["seas"] = Rcpp::wrap(seas_stored);
                // output["seas_accept"] = static_cast<double>(seas_accept / ntotal);

                output["infer_rho"] = rho_prior.infer;
                output["rho"] = Rcpp::wrap(rho_stored);
                // output["rho_accept"] = rho_accept / ntotal;

                output["infer_par1"] = par1_prior.infer;
                output["par1"] = Rcpp::wrap(par1_stored);

                output["infer_par2"] = par2_prior.infer;
                output["par2"] = Rcpp::wrap(par2_stored);
                // output["lag_accept"] = lag_accept / ntotal;

                output["infer_zintercept"] = zintercept_infer;
                if (zintercept_infer)
                {
                    output["zintercept"] = Rcpp::wrap(zintercept_stored);
                }

                output["infer_zzcoef"] = zzcoef_infer;
                if (zzcoef_infer)
                {
                    output["zzcoef"] = Rcpp::wrap(zzcoef_stored);
                }

                output["hmc_accept"] = hmc_accept / ntotal;

                // HMC diagnostics
                Rcpp::List hmc_list = Rcpp::List::create(
                    Rcpp::Named("acceptance_rate") = hmc_accept / static_cast<double>(ntotal),
                    Rcpp::Named("leapfrog_step_size") = hmc_step_size,
                    Rcpp::Named("n_leapfrog") = L,
                    Rcpp::Named("mass_diag") = Rcpp::wrap(mass_diag)
                );
                if (hmc_diag.energy_diff.n_elem > 0)
                {
                    hmc_list["diagnostics"] = hmc_diag.to_list();
                }
                output["hmc"] = hmc_list;
            }

            return output;
        }


        void infer(Model &model, const arma::vec &y, const bool &verbose = VERBOSE)
        {
            const unsigned int nT = y.n_elem - 1;
            model.seas.X = Season::setX(nT, model.seas.period, model.seas.P);
            arma::vec ztmp(nT + 1, arma::fill::ones);
            if (model.zero.inflated)
            {
                z_stored.set_size(nT + 1, nsample);
                z_stored.ones();
                ztmp.zeros();
                for (unsigned int t = 1; t < y.n_elem; t++)
                {
                    if (y.at(t) > EPS)
                    {
                        ztmp.at(t) = 1.;
                    }
                }
            }
            model.zero.setZ(ztmp, nT);

            std::map<std::string, AVAIL::Dist> lag_list = LagDist::lag_list;
            if (!model.dlag.truncated && lag_list[model.dlag.name] == AVAIL::Dist::nbinomp)
            {
                // iterative transfer function
                model.dlag.nL = nT;
                model.dlag.Fphi = LagDist::get_Fphi(model.dlag);
                model.dlag.truncated = true;
            }

            std::map<std::string, AVAIL::Dist> obs_list = ObsDist::obs_list;


            wt = arma::randn(nT + 1) * 0.01;
            wt.at(0) = 0.;
            wt.subvec(1, model.nP) = arma::abs(wt.subvec(1, model.nP));
            psi = arma::cumsum(wt);
            wt_accept.set_size(nT + 1);
            wt_accept.zeros();

            #ifdef DGTF_DO_BOUND_CHECK
            bound_check(wt, "Disturbance::init");
            #endif

            psi_stored.set_size(nT + 1, nsample);
            psi_stored.zeros();
            ApproxDisturbance approx_dlm(nT, model.fgain);

            // ================================================================
            // Initialize HMC adaptation machinery
            // ================================================================
            DualAveraging_arma da_adapter;
            MassAdapter_arma mass_adapter;
            if (update_static)
            {
                // Initialize dual averaging
                da_adapter = DualAveraging_arma(hmc_step_size_init, hmc_target_accept);

                // Initialize diagnostics
                hmc_diag = HMCDiagnostics_arma(ntotal, nburnin);

                // Initialize mass adapter with current unconstrained parameters
                arma::vec params0 = Static::init_eta(param_selected, model, true);
                arma::vec q0 = Static::eta2tilde(
                    params0, param_selected, W_prior.name, par1_prior.name,
                    model.dobs.name, model.seas.period, model.seas.in_state);
                mass_adapter.init(q0);
            }

            for (unsigned int b = 0; b < ntotal; b++)
            {
                Rcpp::checkUserInterrupt();

                if (obs_list[model.dobs.name] == AVAIL::Dist::nbinomp)
                {
                    arma::vec psi_old = psi;
                    psi = Posterior::update_pg_psi(model, y, psi_old);
                    wt.subvec(1, nT) = arma::diff(psi);
                }
                else
                {
                    approx_dlm.set_Fphi(model.dlag, model.dlag.nL);
                    Posterior::update_wt(wt, wt_accept, approx_dlm, y, model, mh_sd);
                    psi = arma::cumsum(wt);
                }

                if (model.zero.inflated)
                {
                    Posterior::update_zt(model, y, wt);
                }


                if (update_static)
                {
                    double energy_diff = 0.0;
                    double grad_norm = 0.0;
                    double accept_prob = Posterior::update_static_hmc(
                        energy_diff, grad_norm,
                        model, y, psi, wt, param_selected,
                        W_prior, seas_prior, rho_prior, par1_prior, par2_prior,
                        hmc_step_size, mass_diag,
                        zintercept_infer, zzcoef_infer, L);

                    hmc_accept += accept_prob;

                    // Store diagnostics
                    if (hmc_diag.energy_diff.n_elem > 0)
                    {
                        hmc_diag.energy_diff.at(b) = energy_diff;
                        hmc_diag.grad_norm.at(b) = grad_norm;
                    }

                    // ========================================================
                    // Dual averaging: adapt step size during burnin
                    // ========================================================
                    if (hmc_dual_averaging && b <= nburnin)
                    {
                        if (b < nburnin)
                        {
                            hmc_step_size = da_adapter.update_step_size(accept_prob);
                        }
                        else if (b == nburnin)
                        {
                            da_adapter.finalize(hmc_step_size, L, hmc_T_target);
                        }

                        // Store adaptation trace
                        if (hmc_diag.step_size_stored.n_elem > 0 && b <= nburnin)
                        {
                            hmc_diag.step_size_stored.at(b) = hmc_step_size;
                            hmc_diag.nleapfrog_stored.at(b) = static_cast<double>(L);
                        }
                    } // if dual averaging

                    // ========================================================
                    // Mass matrix adaptation during burnin
                    // ========================================================
                    if (b < nburnin)
                    {
                        // Collect unconstrained parameter samples for variance estimation
                        arma::vec params_curr = Static::init_eta(param_selected, model, true);
                        arma::vec q_curr = Static::eta2tilde(
                            params_curr, param_selected, W_prior.name, par1_prior.name,
                            model.dobs.name, model.seas.period, model.seas.in_state);
                        mass_adapter.update(q_curr);

                        // Update mass matrix once at midpoint of burnin
                        if (b == nburnin / 2)
                        {
                            mass_diag = mass_adapter.get_mass_diag();

                            // CRITICAL: Reset dual averaging for the new geometry
                            if (hmc_dual_averaging)
                            {
                                da_adapter = DualAveraging_arma(hmc_step_size, hmc_target_accept);
                            }
                        }
                    } // mass matrix adaptation
                } // if update_static


                bool saveiter = b > nburnin && ((b - nburnin - 1) % nthin == 0);
                if (saveiter || b == (ntotal - 1))
                {
                    unsigned int idx_run;
                    if (saveiter)
                    {
                        idx_run = (b - nburnin - 1) / nthin;
                    }
                    else
                    {
                        idx_run = nsample - 1;
                    }

                    // log_marg_stored.at(idx_run) = log_marg;
                    psi_stored.col(idx_run) = psi;
                    if (model.zero.inflated)
                    {
                        z_stored.col(idx_run) = model.zero.z;
                    }

                    if (update_static)
                    {
                        W_stored.at(idx_run) = model.derr.par1;
                        seas_stored.col(idx_run) = model.seas.val;
                        rho_stored.at(idx_run) = model.dobs.par2;
                        par1_stored.at(idx_run) = model.dlag.par1;
                        par2_stored.at(idx_run) = model.dlag.par2;
                        zintercept_stored.at(idx_run) = model.zero.intercept;
                        zzcoef_stored.at(idx_run) = model.zero.coef;
                    }
                }

                if (verbose)
                {
                    Rprintf("\rProgress: %u/%u", b, ntotal - 1);
                }

            } // end a single iteration

            if (verbose)
            {
                Rprintf("\n");
            }

            return;
        }

    private:
        arma::vec log_marg_stored;

        unsigned int L = 50;
        double hmc_step_size = 0.005;
        double hmc_step_size_init = 0.005;

        double hmc_T_target = 2.0;
        double hmc_target_accept = 0.75;
        bool hmc_dual_averaging = true;
        arma::vec mass_diag;                // nparam x 1, diagonal mass matrix
        HMCDiagnostics_arma hmc_diag;       // diagnostics storage

        double mh_sd = 1.0;
        unsigned int nburnin = 100;
        unsigned int nthin = 1;
        unsigned int nsample = 100;
        unsigned int ntotal = 200;
        unsigned int max_lag = 50;

        bool update_static = true;
        double hmc_accept = 0.;
        unsigned int nparam = 1; // number of unknown static parameters
        std::vector<std::string> param_selected = {"W"};

        arma::mat z_stored; // (nT + 1) x nsample
        bool zintercept_infer = false;
        bool zzcoef_infer = false;
        arma::vec zintercept_stored;
        arma::vec zzcoef_stored;

        arma::vec wt, psi;
        arma::vec wt_accept; // nsample x 1
        arma::mat psi_stored; // (nT + 1) x nsample

        Prior seas_prior;
        arma::mat seas_stored; // period x nsample
        // double seas_accept = 0.;

        Prior rho_prior;
        arma::vec rho_stored;
        // double rho_accept = 0.;

        Prior par1_prior;
        arma::vec par1_stored;
        // double lag_accept = 0.;

        Prior par2_prior;
        arma::vec par2_stored;

        Prior W_prior;
        arma::vec W_stored;
        // double W_accept = 0.;
    };
}

#endif
