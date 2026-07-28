#include <chrono>
#include "Model.h"
#include "LinearBayes.h"
#include "SequentialMonteCarlo.h"
#include "MCMC.h"
#include "VariationalBayes.h"
#include "Diagnostics.h"

#include <progress.hpp>
#include <progress_bar.hpp>

using namespace Rcpp;
// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(RcppArmadillo, BH, RcppProgress, pg)]]


// [[Rcpp::export]]
double ess_cpp(const arma::vec& x) {
    return diagnostics::effective_sample_size(x);
}

// [[Rcpp::export]]
double split_rhat_cpp(const arma::mat& X) {
    return diagnostics::split_rhat(X);
}

//' @export
// [[Rcpp::export]]
Rcpp::NumericVector lognorm2serial(const double &mu, const double &sd2)
{
    return lognorm::lognorm2serial(mu, sd2);
}

//' @export
// [[Rcpp::export]]
Rcpp::NumericVector serial2lognorm(const double &m, const double &s2)
{
    return lognorm::serial2lognorm(m, s2);
}


//' @export
// [[Rcpp::export]]
Rcpp::NumericVector dnbinom0(
        const unsigned int &nL,
        const double &kappa,
        const double &r)
{
    return Rcpp::wrap(nbinom::dnbinom(nL, kappa, r));
}


//' @export
// [[Rcpp::export]]
Rcpp::NumericVector dlognorm0(
        const unsigned int &nL,
        const double &mu,
        const double &sd2)
{
    return Rcpp::wrap(lognorm::dlognorm(nL, mu, sd2));
}




//' @export
// [[Rcpp::export]]
Rcpp::List dgtf_default_algo_settings(const std::string &method)
{
    std::map<std::string, AVAIL::Algo> algo_list = AVAIL::algo_list;
    std::string method_name = tolower(method);

    Rcpp::List opts = SMC::SequentialMonteCarlo::default_settings();
    switch (algo_list[method_name])
    {
    case AVAIL::Algo::LinearBayes:
    {
        opts = LBA::LinearBayes::default_settings();
        break;
    }
    case AVAIL::Algo::MCS:
    {
        opts = SMC::MCS::default_settings();
        break;
    }
    case AVAIL::Algo::FFBS:
    {
        opts = SMC::FFBS::default_settings();
        break;
    }
    case AVAIL::Algo::TFS:
    {
        opts = SMC::TFS::default_settings();
        break;
    }
    case AVAIL::Algo::ParticleLearning:
    {
        opts = SMC::PL::default_settings();
        break;
    }
    case AVAIL::Algo::MCMC:
    {
        opts = MCMC::Disturbance::default_settings();
        break;
    }
    case AVAIL::Algo::HybridVariation:
    {
        opts = VB::Hybrid::default_settings();
        break;
    }
    default:
        throw std::invalid_argument("dgtf_default_algo_settings: Invalid algorithm.");
        break;
    }

    return opts;
}

//' @export
// [[Rcpp::export]]
Rcpp::List dgtf_default_model()
{
    return Model::default_settings();
}


//' @export
// [[Rcpp::export]]
Rcpp::List dgtf_simulate(
    const Rcpp::List &settings,
    const unsigned int &ntime,
    const double &y0 = 0.,
    const Rcpp::Nullable<Rcpp::NumericVector> &z = R_NilValue) // p x (<=ntime), zero inflation
{
    std::map<std::string, SysEq::Evolution> sys_list = SysEq::sys_list;
    Model model(settings);
    if (z.isNotNull())
    {
        arma::vec zvec = Rcpp::as<arma::vec>(z);
        model.zero.setZ(zvec, ntime);
    }


    arma::vec theta0(model.nP, arma::fill::zeros);
    if (sys_list[model.fsys] == SysEq::Evolution::identity)
    {
        Rcpp::List param_settings = settings["param"];
        if (!param_settings.containsElementNamed("lag"))
        {
            throw std::invalid_argument("Model::init - autoregressive coefficients are missing.");
        }

        theta0 = Rcpp::as<Rcpp::NumericVector>(param_settings["lag"]); // ar coefficients
    }

    arma::vec psi, ft, lambda, y;
    arma::mat Theta;
    StateSpace::simulate(y, lambda, ft, Theta, psi, model, ntime, y0, theta0, model.derr.full_rank);

    Rcpp::List output;
    output["y"] = Rcpp::wrap(y.t());
    output["nlag"] = model.dlag.nL;
    output["psi"] = Rcpp::wrap(psi.t());
    output["Theta"] = Rcpp::wrap(Theta);
    output["ft"] = Rcpp::wrap(ft.t());
    output["lambda"] = Rcpp::wrap(lambda.t());

    if (model.zero.inflated)
    {
        output["prob"] = Rcpp::wrap(model.zero.prob.t());
        output["z"] = Rcpp::wrap(model.zero.z.t());
    }

    return output;
}


//' @export
// [[Rcpp::export]]
Rcpp::List dgtf_infer(
    const Rcpp::List &model_settings,
    const arma::vec &y_in,
    const std::string &method,
    const Rcpp::List &method_settings,
    const Rcpp::Nullable<Rcpp::NumericMatrix> &X = R_NilValue, // p x (ntime + 1), zero inflation
    const bool &add_y0 = false)
{
    // Model model = dgtf_initialize(model_settings);
    Model model(model_settings);
    std::map<std::string, AVAIL::Algo> algo_list = AVAIL::algo_list;
    std::string algo_name = tolower(method);

    arma::vec y;
    if (add_y0)
    {
        y.set_size(y_in.n_elem + 1);
        y.at(0) = 0.;
        y.tail(y_in.n_elem) = y_in;
    }
    else
    {
        y = y_in;
    }

    const unsigned int nT = y.n_elem - 1;
    model.seas.X = Season::setX(nT, model.seas.period, model.seas.P);

    if (X.isNotNull())
    {
        arma::mat Xzero = Rcpp::as<arma::mat>(X);
        model.zero.setX(Xzero);
    }

    Rcpp::List output;

    switch (algo_list[algo_name])
    {
    case AVAIL::Algo::LinearBayes:
    {
        y.clamp(0.01 / static_cast<double>(model.nP), y.max());
        LBA::LinearBayes linear_bayes(method_settings);
        linear_bayes.filter(model, y);
        linear_bayes.smoother(model, y);
        output = linear_bayes.get_output(model);
        break;
    } // case Linear Bayes
    case AVAIL::Algo::MCS:
    {
        SMC::MCS mcs(model, method_settings);
        mcs.infer(model, y);
        output = mcs.get_output();
        break;
    } // case MCS
    case AVAIL::Algo::FFBS:
    {
        SMC::FFBS ffbs(model, method_settings);
        ffbs.infer(model, y);
        output = ffbs.get_output();
        break;
    } // case FFBS
    case AVAIL::Algo::TFS:
    {
        SMC::TFS tfs(model, method_settings);
        tfs.infer(model, y);
        output = tfs.get_output();
        break;
    } // case TFS
    case AVAIL::Algo::ParticleLearning:
    {
        SMC::PL pl(model, method_settings);
        pl.infer(model, y);
        output = pl.get_output();
        break;
    } // case particle learning
    case AVAIL::Algo::MCMC:
    {
        MCMC::Disturbance mcmc(model, method_settings);
        mcmc.infer(model, y);
        output = mcmc.get_output();
        break;
    }
    case AVAIL::Algo::HybridVariation:
    {
        VB::Hybrid hvb(model, method_settings);
        hvb.infer(model, y);
        output = hvb.get_output();
        break;
    }
    default:
    {
        throw std::invalid_argument("Unknown algorithm " + algo_name);
        break;
    }
    } // switch by algorithm

    Rcpp::List out;
    out["fit"] = output;

    return out;
}




//' @export
// [[Rcpp::export]]
Rcpp::List dgtf_posterior_predictive(
    const Rcpp::List &output, 
    const Rcpp::List &model_opts, 
    const arma::vec &y, 
    const unsigned int &nrep = 0,
    const Rcpp::Nullable<Rcpp::NumericVector> &Rt = R_NilValue,
    const double &level = 0.95
)
{
    std::map<std::string, AVAIL::Dist> obs_list = ObsDist::obs_list;

    const unsigned int ntime = y.n_elem - 1;
    Model model(model_opts);
    model.seas.X = Season::setX(ntime, model.seas.period, model.seas.P);

    arma::mat psi_stored;
    arma::cube Theta_stored;
    bool use_theta = false;
    if (output.containsElementNamed("Theta"))
    {
        // Theta_stored: nP x (nT + 1) x nsample
        Theta_stored = Rcpp::as<arma::cube>(output["Theta"]);
        use_theta = true;
    } 
    else if (output.containsElementNamed("mt"))
    {
        arma::mat mt = Rcpp::as<arma::mat>(output["mt"]); // np x (nT + 1)
        arma::cube Ct = Rcpp::as<arma::cube>(output["Ct"]); // np x np x (nT + 1)

        Theta_stored.set_size(mt.n_rows, mt.n_cols, 3);
        Theta_stored.zeros();
        for (unsigned int t = 0; t < mt.n_cols; t++)
        {
            arma::vec mt_sd = arma::sqrt(arma::vectorise(Ct.slice(t).diag()));
            Theta_stored.slice(0).col(t) = mt.col(t) - 2. * mt_sd;
            Theta_stored.slice(1).col(t) = mt.col(t);
            Theta_stored.slice(2).col(t) = mt.col(t) + 2. * mt_sd;
        }

        use_theta = true;
    }
    else if (output.containsElementNamed("psi_stored"))
    {
        psi_stored = Rcpp::as<arma::mat>(output["psi_stored"]);
    }
    else if (output.containsElementNamed("psi"))
    {
        psi_stored = Rcpp::as<arma::mat>(output["psi"]);
    }
    else if (output.containsElementNamed("psi_filter"))
    {
        psi_stored = Rcpp::as<arma::mat>(output["psi_filter"]);
    }
    else
    {
        throw std::invalid_argument("No state samples.");
    }

    unsigned int nsample;
    if (use_theta)
    {
        nsample = Theta_stored.n_slices;
    }
    else
    {
        nsample = psi_stored.n_cols;
    }

    arma::vec W_stored(nsample);
    W_stored.fill(model.derr.par1);
    if (output.containsElementNamed("W"))
    {
        W_stored = Rcpp::as<arma::vec>(output["W"]);
    }

    arma::vec rho_stored(nsample);
    rho_stored.fill(model.dobs.par2);
    if (output.containsElementNamed("rho"))
    {
        rho_stored = Rcpp::as<arma::vec>(output["rho"]);
    }


    arma::vec par1_stored(nsample);
    par1_stored.fill(model.dlag.par1);
    if (output.containsElementNamed("par1"))
    {
        par1_stored = Rcpp::as<arma::vec>(output["par1"]);
    }


    arma::vec par2_stored(nsample);
    par2_stored.fill(model.dlag.par2);
    if (output.containsElementNamed("par2"))
    {
        par2_stored = Rcpp::as<arma::vec>(output["par2"]);
    }

    arma::mat seas_stored(model.seas.period, nsample, arma::fill::zeros);
    if (output.containsElementNamed("seas"))
    {
        seas_stored = Rcpp::as<arma::mat>(output["seas"]);
    }
    else
    {
        for (unsigned int i = 0; i < nsample; i++)
        {
            seas_stored.col(i) = model.seas.val;
        }
    }

    arma::vec hpsi_true;
    arma::mat hpsi_res(ntime + 1, nsample, arma::fill::zeros);
    arma::mat hpsi_stored(ntime + 1, nsample, arma::fill::zeros);
    if (Rt.isNotNull())
    {
        hpsi_true = Rcpp::as<arma::vec>(Rt);
    }

    arma::cube yhat;
    if (nrep > 0)
    {
        yhat = arma::zeros<arma::cube>(ntime + 1, nsample, nrep);

    }
    arma::vec chi_sqr(nsample, arma::fill::zeros);
    Progress p(nsample*ntime, true);
    #ifdef DGTF_USE_OPENMP
    #pragma omp parallel for num_threads(NUM_THREADS) schedule(runtime)
    #endif
    for (unsigned int i = 0; i < nsample; i++)
    {
        arma::vec hpsi;
        arma::mat Theta;
        if (use_theta)
        {
            Theta = Theta_stored.slice(i); // ar(nP) x (nT + 1)
        }
        else
        {
            arma::vec psi = psi_stored.col(i);
            hpsi = GainFunc::psi2hpsi<arma::vec>(psi, model.fgain);
            hpsi_stored.col(i) = hpsi;
            if (Rt.isNotNull())
            {
                hpsi_res.col(i) = arma::abs(hpsi - hpsi_true);
            }
        }

        Model mod = model;
        mod.dobs.par2 = rho_stored.at(i);
        mod.derr.par1 = W_stored.at(i);
        mod.dlag.par1 = par1_stored.at(i);
        mod.dlag.par2 = par2_stored.at(i);
        mod.seas.val = seas_stored.col(i);

        arma::vec ft(ntime + 1, arma::fill::zeros);
        arma::vec yres2(ntime + 1, arma::fill::zeros);
        arma::vec yvar = yres2;
        for (unsigned int t = 1; t <= ntime; t++)
        {
            double eta;
            if (use_theta)
            {
                ft.at(t) = TransFunc::func_ft(
                    model.ftrans, model.fgain, model.dlag,
                    model.seas, t, Theta.col(t), y);

                eta = ft.at(t);
            }
            else
            {
                ft.at(t) = TransFunc::func_ft(t, y, ft, hpsi, mod.dlag, mod.ftrans);

                eta = ft.at(t);
                if (mod.seas.period > 0)
                {
                    eta += arma::dot(mod.seas.X.col(t), mod.seas.val);
                }
            }

            double lambda = std::abs(LinkFunc::ft2mu(eta, model.flink));
            double mean, var;
            switch (obs_list[model.dobs.name])
            {
            case AVAIL::Dist::nbinomm:
            {
                mean = lambda;
                var = std::abs(nbinomm::var(lambda, model.dobs.par2));
                break;
            }
            case AVAIL::Dist::nbinomp:
            {
                mean = nbinom::mean(lambda, model.dobs.par2);
                var = nbinom::var(lambda, model.dobs.par2);
                break;
            }
            case AVAIL::Dist::poisson:
            {
                mean = lambda;
                var = lambda;
                break;
            }
            default:
            {
                throw std::invalid_argument("Unknown observation distribution.");
                break;
            }
            }


            yres2.at(t) = 2. * std::log(std::abs(y.at(t) - lambda) + EPS);
            yvar.at(t) = std::log(var + EPS);

            if (nrep > 0)
            {
                for (unsigned int j = 0; j < nrep; j++)
                {
                    yhat.at(t, i, j) = ObsDist::sample(lambda, mod.dobs.par2, mod.dobs.name);
                }
            }

            p.increment(); 
        } // End iterations over time

        chi_sqr.at(i) = arma::mean(arma::exp(yres2 - yvar));
    } // End iterations over posterior samples


    Rcpp::List output2;
    output2["chi"] = arma::mean(chi_sqr);

    const double alpha = 1.0 - level;
    arma::vec prob = {alpha / 2.0, 0.5, 1.0 - alpha / 2.0};
    if (nrep > 0)
    {
        arma::cube ytmp = yhat.reshape(ntime + 1, nsample * nrep, 1);
        arma::mat yhat2 = ytmp.slice(0);
        output2["crps"] = calculate_crps(y, yhat2);
        
        try
        {
            arma::mat yqt = arma::quantile(yhat2, prob, 1);
            output2["yhat"] = Rcpp::wrap(yqt);

            // Drop the t=0 initial row to match how `chi` and the time loop are averaged.
            arma::vec y_eval = y.subvec(1, ntime);
            arma::vec lo_y = yqt.col(0).subvec(1, ntime);
            arma::vec hi_y = yqt.col(2).subvec(1, ntime);

            arma::uvec yhat_mask = (y_eval >= lo_y - EPS) % (y_eval <= hi_y + EPS);
            output2["coverage_yhat"] = arma::mean(arma::conv_to<arma::vec>::from(yhat_mask));
            output2["width_yhat"] = arma::mean(hi_y - lo_y);

            // Winkler / Gneiting-Raftery interval score (proper scoring rule).
            arma::vec is_score = (hi_y - lo_y) + (2.0 / alpha) * arma::clamp(lo_y - y_eval, 0.0, arma::datum::inf) + (2.0 / alpha) * arma::clamp(y_eval - hi_y, 0.0, arma::datum::inf);
            output2["interval_score_yhat"] = arma::mean(is_score);
        }
        catch(const std::exception& e)
        {
            REprintf("%s\n", e.what());
            output2["yhat"] = Rcpp::wrap(yhat2);
        }
    }

    arma::mat hpsi_qt = arma::quantile(hpsi_stored, prob, 1);
    output2["Rt"] = Rcpp::wrap(hpsi_qt);

    // Always-available width on the latent intensity:
    arma::vec lo_R = hpsi_qt.col(0).subvec(1, ntime);
    arma::vec hi_R = hpsi_qt.col(2).subvec(1, ntime);
    output2["width_Rt"] = arma::mean(hi_R - lo_R);

    if (!use_theta && Rt.isNotNull())
    {
        arma::vec hpsi_eval = hpsi_true.subvec(1, ntime);

        // Aggregate residual matrix the same way as the existing rmse/mae:
        output2["rmse_Rt"] = std::sqrt(arma::mean(arma::mean(arma::pow(hpsi_res, 2))));
        output2["mae_Rt"] = arma::mean(arma::mean(hpsi_res));

        arma::uvec Rt_mask = (hpsi_eval >= lo_R - EPS) % (hpsi_eval <= hi_R + EPS);
        output2["coverage_Rt"] = arma::mean(arma::conv_to<arma::vec>::from(Rt_mask));

        arma::vec is_R = (hi_R - lo_R) + (2.0 / alpha) * arma::clamp(lo_R - hpsi_eval, 0.0, arma::datum::inf) + (2.0 / alpha) * arma::clamp(hpsi_eval - hi_R, 0.0, arma::datum::inf);
        output2["interval_score_Rt"] = arma::mean(is_R);
    }

    output2["level"] = level;
    return output2;
}



//' @export
// [[Rcpp::export]]
Rcpp::List dgtf_forecast(
    const Rcpp::List &output,
    const Rcpp::List &model_settings,
    const arma::vec &y,   // (nT + 1) x 1
    const unsigned int nrep = 100,
    const unsigned int &k = 1, // k-step-ahead forecasting,
    const Rcpp::Nullable<Rcpp::NumericVector> &ypred_true = R_NilValue,
    const bool &only_evaluate_last_one = false,
    const bool &return_all_samples = false
)
{
    std::map<std::string, AVAIL::Dist> obs_list = ObsDist::obs_list;
    const unsigned int ntime = y.n_elem - 1;
    Model model(model_settings);
    model.seas.X = Season::setX(ntime + k, model.seas.period, model.seas.P); 
    // X: period x (ntime + k + 1)

    // Identity-system (AR) fits store per-time-step coefficient
    // vectors in `Theta`, not a scalar latent driver in `psi_stored`.
    // Forecast for these models needs the cube; flag it here and
    // pull it once outside the per-sample loop.
    std::map<std::string, SysEq::Evolution> sys_list_local =
        SysEq::sys_list;
    const bool is_identity =
        (sys_list_local[model.fsys] == SysEq::Evolution::identity);

    arma::cube Theta_stored;
    if (is_identity)
    {
        if (!output.containsElementNamed("Theta"))
        {
            throw std::invalid_argument(
                "Identity-system (AR) forecast requires `Theta` in the "
                "fit output. HVA fits expose this; MCMC fits currently "
                "do not. Forecasting AR models from MCMC fits is not "
                "yet supported.");
        }
        // nP x (ntime + 1) x nsample
        Theta_stored = Rcpp::as<arma::cube>(output["Theta"]);
    }

    arma::mat psi_stored;
    if (output.containsElementNamed("psi_stored"))
    {
        // ntime x nsample
        psi_stored = Rcpp::as<arma::mat>(output["psi_stored"]);
    }
    else if (output.containsElementNamed("psi"))
    {
        // ntime x 3
        psi_stored = Rcpp::as<arma::mat>(output["psi"]);
    }
    else if (output.containsElementNamed("psi_filter"))
    {
        
        psi_stored = Rcpp::as<arma::mat>(output["psi_filter"]);
    }
    else
    {
        throw std::invalid_argument("No state samples.");
    }

    const unsigned int nsample = psi_stored.n_cols;

    arma::vec W_stored(nsample);
    W_stored.fill(model.derr.par1);
    if (output.containsElementNamed("W"))
    {
        W_stored = Rcpp::as<arma::vec>(output["W"]);
    }
    else if (output.containsElementNamed("W_forecast"))
    {
        // Discounted-variational fit: no inferred scalar W, but the
        // engine stored the last in-sample W_t (the discount-derived
        // variance at t = T). Use it as the per-step random-walk
        // variance for the forecast horizon.
        W_stored.fill(Rcpp::as<double>(output["W_forecast"]));
    }

    arma::vec rho_stored(nsample);
    rho_stored.fill(model.dobs.par2);
    if (output.containsElementNamed("rho"))
    {
        rho_stored = Rcpp::as<arma::vec>(output["rho"]);
    }


    arma::vec par1_stored(nsample);
    par1_stored.fill(model.dlag.par1);
    if (output.containsElementNamed("par1"))
    {
        par1_stored = Rcpp::as<arma::vec>(output["par1"]);
    }


    arma::vec par2_stored(nsample);
    par2_stored.fill(model.dlag.par2);
    if (output.containsElementNamed("par2"))
    {
        par2_stored = Rcpp::as<arma::vec>(output["par2"]);
    }

    arma::mat seas_stored(model.seas.period, nsample, arma::fill::zeros);
    if (output.containsElementNamed("seas"))
    {
        seas_stored = Rcpp::as<arma::mat>(output["seas"]);
    }
    else
    {
        for (unsigned int i = 0; i < nsample; i++)
        {
            seas_stored.col(i) = model.seas.val;
        }
    }

    const unsigned int nsample_eff = is_identity
                                         ? Theta_stored.n_slices
                                         : psi_stored.n_cols;
    const unsigned int nP = is_identity ? Theta_stored.n_rows : 0u;

    arma::cube ypred(ntime + k + 1, nrep, nsample_eff, arma::fill::zeros);
    arma::cube ymean = ypred;
    arma::cube yvar = ypred;

    Progress p(nsample_eff, true);
    #ifdef DGTF_USE_OPENMP
    #pragma omp parallel for num_threads(NUM_THREADS) schedule(runtime)
    #endif
    for (unsigned int i = 0; i < nsample_eff; i++)
    {
        // Per-sample setup (shared across both branches): seed ytmp
        // with the in-sample observations and prepare the per-sample
        // model copy with this iteration's static parameters.
        arma::mat ytmp(ntime + 1 + k, nrep, arma::fill::zeros);
        for (unsigned int s = 0; s < nrep; s++)
            ytmp.submat(0, s, ntime, s) = y;

        Model mod = model;
        mod.dobs.par2 = rho_stored.at(i);
        mod.derr.par1 = W_stored.at(i);
        mod.dlag.par1 = par1_stored.at(i);
        mod.dlag.par2 = par2_stored.at(i);
        mod.seas.val = seas_stored.col(i);

        if (is_identity)
        {
            // ----- Identity-system (AR) forecast -----
            //
            // Each lag-l coefficient component evolves as an
            // independent random walk with variance W per step
            // (or stays constant when W = 0). The initial coefficient
            // vector is the posterior sample's last-time-step state
            // taken from the Theta cube.
            //
            // Transfer function with uniform lag (`Fphi == 1`):
            //     f_t = sum_{l=0}^{nP-1} h(theta_l) * y_{t-1-l}
            arma::mat theta_state(nP, nrep);
            const arma::vec theta_last =
                Theta_stored.slice(i).col(ntime); // nP x 1
            for (unsigned int s = 0; s < nrep; s++)
                theta_state.col(s) = theta_last;

            const double sd_W =
                std::sqrt(std::max(0.0, mod.derr.par1));

            for (unsigned int j = 1; j <= k; j++)
            {
                const unsigned int tidx = ntime + j;
                for (unsigned int s = 0; s < nrep; s++)
                {
                    // Identity dynamics: theta_t = theta_{t-1} + N(0, W * I).
                    if (sd_W > 0.0)
                    {
                        for (unsigned int l = 0; l < nP; l++)
                            theta_state.at(l, s) +=
                                sd_W * arma::randn();
                    }

                    // Transfer function (uniform lag, no Fphi weighting).
                    double f_t = 0.0;
                    const unsigned int nelem = std::min(tidx, nP);
                    for (unsigned int l = 0; l < nelem; l++)
                    {
                        const double psi_l = theta_state.at(l, s);
                        const double hpsi_l =
                            GainFunc::psi2hpsi(psi_l, model.fgain);
                        const double y_lag = ytmp.at(tidx - 1 - l, s);
                        f_t += hpsi_l * y_lag;
                    }

                    double eta = f_t;
                    if (mod.seas.period > 0)
                        eta += arma::dot(mod.seas.X.col(tidx),
                                         mod.seas.val);

                    const double lambda =
                        LinkFunc::ft2mu(eta, model.flink);
                    double mean = 0.0, var = 0.0;
                    switch (obs_list[model.dobs.name])
                    {
                    case AVAIL::Dist::nbinomm:
                        mean = lambda;
                        var = std::abs(
                            nbinomm::var(lambda, model.dobs.par2));
                        break;
                    case AVAIL::Dist::nbinomp:
                        mean = nbinom::mean(lambda, model.dobs.par2);
                        var = nbinom::var(lambda, model.dobs.par2);
                        break;
                    case AVAIL::Dist::poisson:
                        mean = lambda;
                        var = lambda;
                        break;
                    default:
                        throw std::invalid_argument(
                            "Unknown observation distribution.");
                    }

                    ymean.at(tidx, s, i) = mean;
                    yvar.at(tidx, s, i) = var;
                    ypred.at(tidx, s, i) = ObsDist::sample(
                        lambda, mod.dobs.par2, mod.dobs.name);
                    ytmp.at(tidx, s) = ypred.at(tidx, s, i);
                }
            }
        }
        else
        {
            // ----- Shift-system forecast (existing path) -----
            //
            // psi_stored is (ntime + 1) x nsample; column i is the
            // per-sample latent-driver trajectory. Continue forward
            // by holding the last value (matches the original behavior).
            arma::mat ft = ytmp;
            arma::mat psi = ytmp;
            const double sd_W = std::sqrt(std::max(0.0, mod.derr.par1));
            for (unsigned int s = 0; s < nrep; s++)
            {
                psi.submat(0, s, ntime, s) = psi_stored.col(i);
                // psi.submat(ntime + 1, s, ntime + k, s)
                //     .fill(psi.at(ntime, s));
                for (unsigned int j = 1; j <= k; j++)
                {
                    psi.at(ntime + j, s) = psi.at(ntime + j - 1, s) + (sd_W > 0.0 ? sd_W * arma::randn() : 0.0);
                }
            }

            arma::mat hpsi =
                GainFunc::psi2hpsi<arma::mat>(psi, model.fgain);

            for (unsigned int j = 1; j <= k; j++)
            {
                const unsigned int tidx = ntime + j;
                for (unsigned int s = 0; s < nrep; s++)
                {
                    ft.at(tidx, s) = TransFunc::func_ft(
                        tidx, ytmp.col(s), ft.col(s), hpsi.col(s),
                        mod.dlag, mod.ftrans);

                    double eta = ft.at(tidx, s);
                    if (mod.seas.period > 0)
                        eta += arma::dot(mod.seas.X.col(tidx),
                                         mod.seas.val);

                    const double lambda =
                        LinkFunc::ft2mu(eta, model.flink);
                    double mean = 0.0, var = 0.0;
                    switch (obs_list[model.dobs.name])
                    {
                    case AVAIL::Dist::nbinomm:
                        mean = lambda;
                        var = std::abs(
                            nbinomm::var(lambda, model.dobs.par2));
                        break;
                    case AVAIL::Dist::nbinomp:
                        mean = nbinom::mean(lambda, model.dobs.par2);
                        var = nbinom::var(lambda, model.dobs.par2);
                        break;
                    case AVAIL::Dist::poisson:
                        mean = lambda;
                        var = lambda;
                        break;
                    default:
                        throw std::invalid_argument(
                            "Unknown observation distribution.");
                    }

                    ymean.at(tidx, s, i) = mean;
                    yvar.at(tidx, s, i) = var;
                    ypred.at(tidx, s, i) = ObsDist::sample(
                        lambda, mod.dobs.par2, mod.dobs.name);
                    ytmp.at(tidx, s) = ypred.at(tidx, s, i);
                }
            }
        }

        p.increment();
    }

    arma::mat ypred_mat(
        ypred.memptr(),
        ypred.n_rows,
        ypred.n_cols * ypred.n_slices,
        /*copy_aux_mem=*/false,
        /*strict=*/true);

    arma::mat ypred_sub;
    if (only_evaluate_last_one)
    {
        ypred_sub = ypred_mat.tail_rows(1);
    }
    else
    {
        ypred_sub = ypred_mat.tail_rows(k);
    }

    arma::vec prob = {0.025, 0.5, 0.975};
    arma::mat yqt = arma::quantile(ypred_sub, prob, 1);

    Rcpp::List out;
    out["yqt"] = Rcpp::wrap(yqt);

    if (return_all_samples)
    {
        out["ypred_samples"] = Rcpp::wrap(ypred_sub);
    }


    if (!ypred_true.isNull())
    {
        arma::mat ymean_mat(
            ymean.memptr(),
            ymean.n_rows,
            ymean.n_cols * ymean.n_slices,
            /*copy_aux_mem=*/false,
            /*strict=*/true);

        arma::mat yvar_mat(
            yvar.memptr(),
            yvar.n_rows,
            yvar.n_cols * yvar.n_slices,
            /*copy_aux_mem=*/false,
            /*strict=*/true);

        arma::vec ypredt;
        if (only_evaluate_last_one)
        {
            ypredt = Rcpp::as<arma::vec>(ypred_true).tail(1);
        }
        else
        {
            ypredt = Rcpp::as<arma::vec>(ypred_true).tail(k);
        }

        out["ypred_true"] = Rcpp::wrap(ypredt);

        arma::uvec ypred_mask = (ypredt >= (yqt.col(0) - EPS)) % (ypredt <= (yqt.col(2) + EPS));
        out["coverage"] = static_cast<double>(arma::accu(ypred_mask)) / static_cast<double>(ypred_mask.n_elem);

        arma::mat ymean_sub;
        arma::mat yvar_sub;
        if (only_evaluate_last_one)
        {
            ymean_sub = ymean_mat.tail_rows(1);
            yvar_sub = yvar_mat.tail_rows(1);
        }
        else
        {
            ymean_sub = ymean_mat.tail_rows(k);
            yvar_sub = yvar_mat.tail_rows(k);
        }
        arma::mat log_diff = arma::log(arma::abs(ymean_sub.each_col() - ypredt) + EPS);
        arma::mat log_var = arma::log(yvar_sub + EPS);
        double chi_sqr = arma::mean(arma::mean(arma::exp(log_diff - log_var)));
        out["chi"] = chi_sqr;

        double crps = calculate_crps(ypredt, ypred_sub);
        out["crps"] = crps;
    }

    return out;
}
