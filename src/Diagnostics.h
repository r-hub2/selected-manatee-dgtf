#pragma once
#include <RcppArmadillo.h>

namespace diagnostics {

inline double effective_sample_size(const arma::vec& draws_in)
{
    arma::vec vals = draws_in.elem(arma::find_finite(draws_in));
    const arma::uword n = vals.n_elem;
    if (n < 2) return static_cast<double>(n);

    vals -= arma::mean(vals);

    arma::uword nfft = 1;
    while (nfft < 2 * n) nfft <<= 1;

    arma::vec padded(nfft, arma::fill::zeros);
    padded.head(n) = vals;

    arma::cx_vec freq = arma::fft(padded);
    freq %= arma::conj(freq);                  // |X(ω)|²
    arma::cx_vec acov_c = arma::ifft(freq);

    arma::vec acov(n);
    for (arma::uword k = 0; k < n; ++k)
        acov(k) = acov_c(k).real() / static_cast<double>(n - k);  // unbiased

    const double gamma0 = acov(0);
    if (!std::isfinite(gamma0) || gamma0 <= 0.0)
        return static_cast<double>(n);

    double sum_rho = 0.0;
    for (arma::uword lag = 1; lag + 1 < n; lag += 2) {
        double pair = acov(lag) / gamma0 + acov(lag + 1) / gamma0;
        if (!std::isfinite(pair) || pair < 0.0) break;             // Geyer IPS
        sum_rho += pair;
    }

    double ess = static_cast<double>(n) / (1.0 + 2.0 * sum_rho);
    return std::min(static_cast<double>(n), std::max(1.0, ess));
} // effective_sample_size

inline double split_rhat(const arma::mat& draws_in)
{
    const arma::uword n_in = draws_in.n_rows;
    const arma::uword m_in = draws_in.n_cols;
    if (m_in < 1)
        Rcpp::stop("split_rhat: at least one chain (column) is required.");
    if (n_in < 4)
        Rcpp::stop("split_rhat: chain length < 4 leaves nothing to split.");

    // Sequentially split each chain in half (drop the middle row when n is odd).
    const arma::uword n_half = n_in / 2;
    const arma::uword m      = 2 * m_in;
    arma::mat X(n_half, m);
    for (arma::uword j = 0; j < m_in; ++j) {
        X.col(2 * j)     = draws_in.col(j).head(n_half);
        X.col(2 * j + 1) = draws_in.col(j).tail(n_half);
    }

    // Drop rows with any non-finite entry.
    arma::uvec good(X.n_rows);
    arma::uword n_good = 0;
    for (arma::uword i = 0; i < X.n_rows; ++i)
        if (X.row(i).is_finite())
            good(n_good++) = i;
    if (n_good < 2)
        Rcpp::stop("split_rhat: fewer than 2 finite rows after filtering.");
    X = X.rows(good.head(n_good));

    const double n  = static_cast<double>(X.n_rows);
    const double mD = static_cast<double>(m);

    arma::rowvec chain_means = arma::mean(X, 0);          // 1 x m
    arma::rowvec chain_vars(m);
    for (arma::uword j = 0; j < m; ++j) {
        arma::vec c = X.col(j) - chain_means(j);
        chain_vars(j) = arma::dot(c, c) / (n - 1.0);
    }

    const double W      = arma::mean(chain_vars);
    const double mu_bar = arma::mean(chain_means);
    // Stan / BDA3 unbiased between-chain variance: divide by (m-1).
    const double B =
        n * arma::sum(arma::square(chain_means - mu_bar)) / (mD - 1.0);

    if (W <= 0.0)
        return (B <= 0.0) ? 1.0 : std::numeric_limits<double>::infinity();

    const double var_hat = ((n - 1.0) / n) * W + (1.0 / n) * B;
    return std::sqrt(var_hat / W);
} // split_rhat

} // namespace diagnostics
