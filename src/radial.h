#ifndef BEHAVSSF_RADIAL_H
#define BEHAVSSF_RADIAL_H

// The order matches .step_families in R/families.R. Derivatives are on
// the reported parameter scale. The objective applies the working-scale link.
enum StepFamily { canonical, gamma_step, lognormal_step, inverse_gaussian_step,
                  halfnormal_step, weibull_step, exponential_step };

struct RadialDensity { double log_density, first_score, second_score; };

inline RadialDensity radial_density(int family, double r, double u,
                                    const double* beta, double gamma_digamma) {
  const double a = beta[0], log_2pi = std::log(2 * std::acos(-1.0));
  switch (family) {
  case gamma_step: {
    const double rate = beta[1];
    return {Rf_dgamma(r, a, 1 / rate, 1),
      std::log(rate) + u - gamma_digamma, a / rate - r};
  }
  case lognormal_step: {
    const double sd = beta[1], z = (u - a) / sd;
    return {-u - std::log(sd) - log_2pi / 2 - z * z / 2, z / sd, (z * z - 1) / sd};
  }
  case inverse_gaussian_step: {
    const double shape = beta[1], delta = std::expm1(u - std::log(a));
    const double scaled = delta / std::sqrt(r), square = scaled * scaled;
    return {(std::log(shape) - log_2pi - 3 * u - shape * square) / 2,
      shape * delta / (a * a), (1 / shape - square) / 2};
  }
  case halfnormal_step: {
    const double z = r / a;
    return {std::log(2.0) - log_2pi / 2 - std::log(a) - z * z / 2,
      (z * z - 1) / a, 0};
  }
  case weibull_step: {
    const double scale = beta[1], log_ratio = u - std::log(scale);
    const double exponent = a * log_ratio, power = std::exp(exponent);
    return {std::log(a) - u + exponent - power,
      1 / a + log_ratio * (1 - power), a * (power - 1) / scale};
  }
  case exponential_step:
    return {std::log(a) - a * r, 1 / a - r, 0};
  default:
    throw std::runtime_error("Unknown step-length family.");
  }
}

inline double radial_draw(int family, double a, double b) {
  switch (family) {
  case gamma_step:
    return Rf_rgamma(a, 1 / b);
  case lognormal_step:
    return Rf_rlnorm(a, b);
  case inverse_gaussian_step: {
    // Michael-Schucany-Haas draw. Rationalising the smaller root avoids
    // cancellation when the squared normal variate is large.
    const double z = norm_rand();
    const double y = (a / b) * (z * z / 2);
    const double ratio = 1 + y + std::sqrt(y) * std::sqrt(2 + y);
    return unif_rand() < 1 / (1 + 1 / ratio) ? a / ratio : a * ratio;
  }
  case halfnormal_step:
    return std::abs(Rf_rnorm(0, a));
  case weibull_step:
    return Rf_rweibull(a, b);
  case exponential_step:
    return Rf_rexp(1 / a);
  default:
    throw std::runtime_error("Unknown step-length family.");
  }
}

#endif
