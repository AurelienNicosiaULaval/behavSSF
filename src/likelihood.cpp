#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Utils.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <memory>
#include <vector>
#include "diel.h"

namespace {
double log_sum(const std::vector<double>& x) {
  double m = *std::max_element(x.begin(), x.end());
  if (!std::isfinite(m)) throw std::runtime_error("Nonfinite log sum.");
  double s = 0;
  for (double v : x) s += std::exp(v - m);
  return m + std::log(s);
}

// Small dense solve with partial pivoting; matrices use row-major storage.
std::vector<double> solve(std::vector<double> a, std::vector<double> b, int k) {
  for (int i = 0; i < k; ++i) {
    int pivot = i;
    for (int j = i + 1; j < k; ++j)
      if (std::abs(a[j*k+i]) > std::abs(a[pivot*k+i])) pivot = j;
    if (std::abs(a[pivot*k+i]) < 1e-14)
      throw std::runtime_error("Stationary system is singular or ill conditioned.");
    for (int j = i; j < k; ++j) std::swap(a[i*k+j], a[pivot*k+j]);
    std::swap(b[i], b[pivot]);
    for (int j = i + 1; j < k; ++j) {
      double f = a[j*k+i] / a[i*k+i];
      for (int l = i; l < k; ++l) a[j*k+l] -= f * a[i*k+l];
      b[j] -= f * b[i];
    }
  }
  std::vector<double> answer(k);
  for (int i = k - 1; i >= 0; --i) {
    double v = b[i];
    for (int j = i + 1; j < k; ++j) v -= a[i*k+j] * answer[j];
    answer[i] = v / a[i*k+i];
  }
  return answer;
}
}

SEXP objective_impl(SEXP theta_, SEXP x_, SEXP logw_,
                    SEXP starts_, SEXP sizes_, SEXP segments_, SEXP k_, SEXP exact_) {
  try {
    const int k = Rf_asInteger(k_);
    const int nrow = INTEGER(Rf_getAttrib(x_, R_DimSymbol))[0];
    const int p = INTEGER(Rf_getAttrib(x_, R_DimSymbol))[1];
    const bool exact = exact_ != R_NilValue;
    const int n = exact ? nrow : Rf_length(starts_);
    const bool diel = exact && Rf_length(exact_) >= 6;
    const int q = diel ? Rf_length(VECTOR_ELT(exact_, 5))/24 : 1;
    const bool hourly_length = exact && Rf_length(exact_) == 7;
    const int r = hourly_length ? Rf_length(VECTOR_ELT(exact_, 6))/24 : 0;
    const int step_start = p*k + k*(k-1)*q;
    const int npar = step_start + k*r;
    if (k < 1 || k > 6 || Rf_length(theta_) != npar || n < 1)
      throw std::runtime_error("Invalid internal objective dimensions.");
    const double* working = REAL(theta_);
    std::vector<double> natural(working, working + Rf_length(theta_));
    std::vector<double> jacobian(Rf_length(theta_), 1);
    int step_index = -1, log_index = -1, persistence_index = -1, target_index = -1;
    std::vector<int> ar_columns, delays;
    std::vector<double> ar_cross(k, 0);
    double shape_offset = 0;
    const double* psi = nullptr;
    const double* active = nullptr;
    if (exact) {
      const int* indices = INTEGER(VECTOR_ELT(exact_, 0));
      if (Rf_length(VECTOR_ELT(exact_, 0)) != 10)
        throw std::runtime_error("Invalid exact kernel indices.");
      step_index = indices[0]; log_index = indices[1];
      persistence_index = indices[2]; target_index = indices[3];
      for (int i = 4; i < 10; ++i) if (indices[i] >= 0) {
        if (indices[i] >= p || indices[4] < 0)
          throw std::runtime_error("Invalid exact AR kernel context.");
        ar_columns.push_back(indices[i]);
        delays.push_back(i == 4 ? 1 : 22+i-5);
      }
      shape_offset = Rf_asReal(VECTOR_ELT(exact_, 3));
      psi = REAL(VECTOR_ELT(exact_, 1));
      active = REAL(VECTOR_ELT(exact_, 2));
      if (step_index < 0 || log_index < 0 || persistence_index < 0 ||
          step_index >= p || log_index >= p || persistence_index >= p || target_index >= p ||
          Rf_length(VECTOR_ELT(exact_,1)) != n || Rf_length(VECTOR_ELT(exact_,2)) != n)
        throw std::runtime_error("Invalid exact kernel context.");
      for (int state = 0; state < k; ++state) {
        const int si = step_index+p*state, li = log_index+p*state;
        const double rate = std::exp(working[si]), shape = std::exp(working[li]);
        natural[si] = -rate; natural[li] = shape-shape_offset;
        jacobian[si] = -rate; jacobian[li] = shape;
        if (!(rate > 0) || !(natural[li]+shape_offset > 0) ||
            !std::isfinite(rate) || !std::isfinite(shape))
          throw std::runtime_error("Nonfinite or nonpositive exact gamma parameters.");
        const int n_lags = ar_columns.size();
        if (n_lags > 0) {
          const int ai = ar_columns[0]+p*state;
          if (n_lags == 2) {
            const int di = ar_columns[1]+p*state;
            const double plus = std::tanh(working[ai]+working[di]);
            const double minus = std::tanh(working[ai]-working[di]);
            natural[ai] = (plus+minus)/2;
            natural[di] = (plus-minus)/2;
            if (!(std::abs(natural[ai])+std::abs(natural[di]) < 1))
              throw std::runtime_error("The absolute AR coefficients must sum to less than one.");
            jacobian[ai] = jacobian[di] = (2-plus*plus-minus*minus)/2;
            ar_cross[state] = (minus*minus-plus*plus)/2;
          } else if (n_lags > 2) {
            // Normalize the hyperbolic sine by the sum of coshes to keep
            // the coefficients strictly inside the unit L1 ball.
            double scale = 0;
            for (int column : ar_columns) scale = std::max(scale, std::abs(working[column+p*state]));
            double total = 0, absolute_sum = 0;
            for (int column : ar_columns) {
              const int j = column+p*state;
              const double pos = std::exp(working[j]-scale), neg = std::exp(-working[j]-scale);
              natural[j] = pos-neg;
              jacobian[j] = pos+neg;
              total += jacobian[j];
            }
            for (int column : ar_columns) {
              const int j = column+p*state;
              natural[j] /= total;
              jacobian[j] = jacobian[j]/total-natural[j]*natural[j];
              absolute_sum += std::abs(natural[j]);
            }
            if (!(absolute_sum < 1))
              throw std::runtime_error("The absolute AR coefficients must sum to less than one.");
          } else {
            const double rho = std::tanh(working[ai]);
            if (!(std::abs(rho) < 1))
              throw std::runtime_error("AR coefficient outside its open parameter range.");
            natural[ai] = rho;
            jacobian[ai] = 1-rho*rho;
          }
        }
      }
    }
    const double* theta = natural.data();
    const double* x = REAL(x_);
    const double* logw = exact ? nullptr : REAL(logw_);
    const int* starts = exact ? nullptr : INTEGER(starts_);
    const int* sizes = exact ? nullptr : INTEGER(sizes_);
    const int* segments = INTEGER(segments_);
    std::vector<double> P(k*k), logP(k*k), delta(k, 1), A(k*k);
    int index = p*k;
    std::unique_ptr<DielTransitions> schedule;
    const int* hour = nullptr;
    if (diel) {
      if (k != 2 || q < 1 || Rf_length(VECTOR_ELT(exact_, 4)) != n)
        throw std::runtime_error("Invalid diel transition context.");
      hour = INTEGER(VECTOR_ELT(exact_, 4));
      for (int t = 0; t < n; ++t) if (hour[t] < 0 || hour[t] > 23)
        throw std::runtime_error("Invalid hour in diel transition context.");
      schedule.reset(new DielTransitions(REAL(VECTOR_ELT(exact_, 5)), q, theta+p*k));
    } else {
      for (int i = 0; i < k; ++i) {
        std::vector<double> logits(k, 0);
        for (int j = 0; j < k; ++j) if (j != i) logits[j] = theta[index++];
        double normalizer = log_sum(logits);
        for (int j = 0; j < k; ++j) {
          logP[i*k+j] = logits[j] - normalizer;
          P[i*k+j] = std::exp(logP[i*k+j]);
        }
      }
      if (k > 1) {
        for (int i = 0; i < k; ++i) for (int j = 0; j < k; ++j)
          A[i*k+j] = i == k-1 ? 1 : P[j*k+i] - (i == j ? 1 : 0);
        std::vector<double> b(k, 0); b[k-1] = 1;
        delta = solve(A, b, k);
      }
      for (double d : delta) if (!(d > 0) || !std::isfinite(d))
        throw std::runtime_error("Invalid stationary probabilities.");
    }
    const double* step_basis = hourly_length ? REAL(VECTOR_ELT(exact_, 6)) : nullptr;
    std::vector<double> step_shift(hourly_length ? 24*k : 0, 0);
    std::vector<double> step_score(hourly_length ? n*k : 0);
    if (hourly_length) for (int h = 0; h < 24; ++h) for (int i = 0; i < k; ++i)
      for (int j = 0; j < r; ++j) step_shift[h*k+i] += step_basis[h+24*j]*theta[step_start+i*r+j];
    const bool conditional_rate = !ar_columns.empty() || hourly_length;
    std::vector<bool> ar_column(p, false);
    for (int column : ar_columns) ar_column[column] = true;
    std::vector<double> emission(n*k), derivative(n*k*p), alpha(n*k),
      backward(n*k, 0), scale(n), gradient(npar, 0), initial_score(k, 0);
    std::vector<double> radial_log_z(k), radial_mean(k), radial_log_mean(k),
      radial_log_rate(k), radial_log_anchor(k), base_angular_log_z(k), base_angular_score(k);
    if (exact) for (int state = 0; state < k; ++state) {
      const double rate = -theta[step_index+p*state];
      const double shape = theta[log_index+p*state]+shape_offset;
      radial_log_z[state] = Rf_lgammafn(shape)-shape*std::log(rate);
      radial_mean[state] = shape/rate;
      radial_log_mean[state] = Rf_digamma(shape)-std::log(rate);
      radial_log_rate[state] = std::log(rate);
      radial_log_anchor[state] = std::log(shape)-radial_log_rate[state];
      const double a = theta[persistence_index+p*state], concentration = std::abs(a);
      const double i0 = Rf_bessel_i(concentration,0,2);
      base_angular_log_z[state] = std::log(2*std::acos(-1.0))+std::log(i0)+concentration;
      base_angular_score[state] = concentration < 1e-5 ? a*(.5-concentration*concentration/16) :
        (a/concentration)*Rf_bessel_i(concentration,1,2)/i0;
    }
    int segment_age = 0;
    for (int s = 0; s < n; ++s) {
      if (s % 512 == 0) R_CheckUserInterrupt();
      segment_age = s > 0 && segments[s] == segments[s-1] ? segment_age+1 : 0;
      const int first = exact ? s : starts[s], count = exact ? 1 : sizes[s];
      for (int state = 0; state < k; ++state) {
        if (exact) {
          const double a = theta[persistence_index+p*state];
          const double bias = target_index < 0 ? 0 : theta[target_index+p*state];
          const double b = bias*active[s];
          const double cos_psi = std::cos(psi[s]);
          const double vx = a+b*cos_psi, vy = b*std::sin(psi[s]);
          const double concentration = std::hypot(vx, vy);
          if (!std::isfinite(concentration))
            throw std::runtime_error("Nonfinite exact angular concentration.");
          const bool no_target = target_index < 0 || active[s] == 0;
          const double scaled_i0 = no_target ? 1 : Rf_bessel_i(concentration, 0, 2);
          const double angular_log_z = no_target ? base_angular_log_z[state] :
            std::log(2*std::acos(-1.0))+std::log(scaled_i0)+concentration;
          double eta = 0;
          for (int v = 0; v < p; ++v) {
            if (conditional_rate && (ar_column[v] || v == step_index || v == log_index)) continue;
            eta += x[s+nrow*v]*theta[v+p*state];
            derivative[(s*k+state)*p+v] = x[s+nrow*v];
          }
          if (!conditional_rate) {
            emission[s*k+state] = eta-(radial_log_z[state]+angular_log_z);
            derivative[(s*k+state)*p+step_index] -= radial_mean[state];
            derivative[(s*k+state)*p+log_index] -= radial_log_mean[state];
          } else {
            double shift = 0, rho_sum = 0;
            std::array<double, 6> ratio;
            std::size_t available = 0;
            for (; available < delays.size() && segment_age >= delays[available]; ++available) {
              const int column = ar_columns[available];
              const double rho = theta[column+p*state];
              ratio[available] = x[s+nrow*column]-radial_log_anchor[state];
              shift -= rho*ratio[available];
              rho_sum += rho;
            }
            if (hourly_length) shift -= step_shift[hour[s]*k+state];
            const double beta_step = theta[step_index+p*state];
            const double rate = shift != 0 ? std::exp(radial_log_rate[state]+shift) : -beta_step;
            if (!(rate > 0) || !std::isfinite(rate))
              throw std::runtime_error("Nonfinite or nonpositive conditional gamma rate.");
            const double step = x[s+nrow*step_index];
            const double log_step = x[s+nrow*log_index];
            const double shape = theta[log_index+p*state]+shape_offset;
            const double rate_score = std::fma(-rate, step, shape);
            // Use the deviance-based gamma density to avoid cancellation at large shapes.
            emission[s*k+state] = Rf_dgamma(step, shape, 1/rate, 1) + eta-angular_log_z;
            if (shape_offset == 2) emission[s*k+state] -= log_step;
            derivative[(s*k+state)*p+step_index] = rate_score*(1-rho_sum)/beta_step;
            derivative[(s*k+state)*p+log_index] = log_step-radial_log_mean[state]+shift+rate_score*rho_sum/shape;
            for (std::size_t j = 0; j < available; ++j)
              derivative[(s*k+state)*p+ar_columns[j]] = -rate_score*ratio[j];
            if (hourly_length) step_score[s*k+state] = -rate_score;
          }
          // I1(k)/(k I0(k)) has limit 1/2 at a zero resultant.
          const double ratio = no_target ? 0 : (concentration < 1e-5 ?
            0.5-concentration*concentration/16 :
            Rf_bessel_i(concentration, 1, 2)/(concentration*scaled_i0));
          derivative[(s*k+state)*p+persistence_index] -= no_target ?
            base_angular_score[state] : ratio*vx;
          if (target_index >= 0)
            derivative[(s*k+state)*p+target_index] -= active[s]*ratio*(b+a*cos_psi);
          if (!std::isfinite(emission[s*k+state]))
            throw std::runtime_error("Nonfinite exact emission density.");
          continue;
        }
        std::vector<double> weighted(count);
        double case_eta = 0;
        for (int row = 0; row < count; ++row) {
          double eta = 0;
          for (int v = 0; v < p; ++v) eta += x[first+row+nrow*v] * theta[v+p*state];
          if (row == 0) case_eta = eta;
          weighted[row] = eta - logw[first+row];
        }
        const double denominator = log_sum(weighted);
        emission[s*k+state] = case_eta - denominator + std::log(count-1.0);
        for (int v = 0; v < p; ++v) {
          double mean = 0;
          for (int row = 0; row < count; ++row)
            mean += std::exp(weighted[row]-denominator) * x[first+row+nrow*v];
          derivative[(s*k+state)*p+v] = x[first+nrow*v] - mean;
        }
      }
    }
    double loglik = 0;
    for (int t = 0; t < n; ++t) {
      std::vector<double> next(k), temp(k);
      const bool first = t == 0 || segments[t] != segments[t-1];
      const double* start = diel ? schedule->start_at(hour[t]) : delta.data();
      const double* log_transition = diel ? schedule->log_at(hour[t]) : logP.data();
      for (int j = 0; j < k; ++j) {
        if (first) next[j] = std::log(start[j]) + emission[t*k+j];
        else {
          for (int i = 0; i < k; ++i) temp[i] = alpha[(t-1)*k+i] + log_transition[i*k+j];
          next[j] = log_sum(temp) + emission[t*k+j];
        }
      }
      scale[t] = log_sum(next);
      loglik += scale[t];
      for (int j = 0; j < k; ++j) alpha[t*k+j] = next[j] - scale[t];
    }
    for (int t = n-2; t >= 0; --t) if (segments[t] == segments[t+1]) {
      const double* log_transition = diel ? schedule->log_at(hour[t+1]) : logP.data();
      for (int i = 0; i < k; ++i) {
        std::vector<double> v(k);
        for (int j = 0; j < k; ++j)
          v[j] = log_transition[i*k+j] + emission[(t+1)*k+j] + backward[(t+1)*k+j];
        backward[t*k+i] = log_sum(v) - scale[t+1];
      }
    }
    for (int t = 0; t < n; ++t) {
      std::vector<double> lp(k), posterior(k);
      for (int i = 0; i < k; ++i) lp[i] = alpha[t*k+i] + backward[t*k+i];
      double normalizer = log_sum(lp);
      for (int i = 0; i < k; ++i) {
        posterior[i] = std::exp(lp[i]-normalizer);
        for (int v = 0; v < p; ++v)
          gradient[v+p*i] += posterior[i] * derivative[(t*k+i)*p+v];
        if (hourly_length) for (int j = 0; j < r; ++j)
          gradient[step_start+i*r+j] += posterior[i]*step_score[t*k+i]*step_basis[hour[t]+24*j];
      }
      const bool first = t == 0 || segments[t] != segments[t-1];
      if (first) {
        if (diel) schedule->add_start_score(hour[t], posterior, gradient.data()+p*k);
        else for (int i = 0; i < k; ++i) initial_score[i] += posterior[i] / delta[i];
      } else {
        index = p*k;
        const double* log_transition = diel ? schedule->log_at(hour[t]) : logP.data();
        const double* transition = diel ? schedule->at(hour[t]) : P.data();
        for (int i = 0; i < k; ++i) {
          std::vector<double> xi(k);
          double row_mass = 0;
          for (int j = 0; j < k; ++j) {
            xi[j] = std::exp(alpha[(t-1)*k+i] + log_transition[i*k+j] +
                              emission[t*k+j] + backward[t*k+j] - scale[t]);
            row_mass += xi[j];
          }
          for (int j = 0; j < k; ++j) if (j != i) {
            double residual = xi[j] - row_mass * transition[i*k+j];
            if (diel) schedule->add_transition_score(hour[t], i, residual, gradient.data()+p*k);
            else gradient[index++] += residual;
          }
        }
      }
    }
    // The initial stationary distribution also depends on the transition logits.
    index = p*k;
    if (!diel) for (int i = 0; i < k; ++i) for (int j = 0; j < k; ++j) if (j != i) {
      std::vector<double> rhs(k);
      for (int l = 0; l < k; ++l)
        rhs[l] = -delta[i] * P[i*k+l] * ((l == j ? 1 : 0) - P[i*k+j]);
      rhs[k-1] = 0;
      std::vector<double> dd = solve(A, rhs, k);
      for (int l = 0; l < k; ++l) gradient[index] += initial_score[l] * dd[l];
      ++index;
    }
    SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));
    SEXP names = PROTECT(Rf_allocVector(STRSXP, 2));
    SEXP grad = PROTECT(Rf_allocVector(REALSXP, npar));
    for (int i = 0; i < npar; ++i) REAL(grad)[i] = -gradient[i]*jacobian[i];
    if (ar_columns.size() == 2) for (int state = 0; state < k; ++state) {
      const int ai = ar_columns[0]+p*state, di = ar_columns[1]+p*state;
      REAL(grad)[ai] -= ar_cross[state]*gradient[di];
      REAL(grad)[di] -= ar_cross[state]*gradient[ai];
    }
    if (ar_columns.size() > 2) for (int state = 0; state < k; ++state)
      for (int a : ar_columns) for (int b : ar_columns) if (a != b) {
        const int i = a+p*state, j = b+p*state;
        REAL(grad)[i] += natural[i]*natural[j]*gradient[j];
      }
    SET_VECTOR_ELT(out, 0, Rf_ScalarReal(-loglik));
    SET_VECTOR_ELT(out, 1, grad);
    SET_STRING_ELT(names, 0, Rf_mkChar("value"));
    SET_STRING_ELT(names, 1, Rf_mkChar("gradient"));
    Rf_setAttrib(out, R_NamesSymbol, names);
    UNPROTECT(3);
    return out;
  } catch (const std::exception& e) { Rf_error("%s", e.what()); }
  return R_NilValue;
}

extern "C" SEXP behavSSF_objective(SEXP theta, SEXP x, SEXP logw,
                                    SEXP starts, SEXP sizes, SEXP segments, SEXP k) {
  return objective_impl(theta,x,logw,starts,sizes,segments,k,R_NilValue);
}
extern "C" SEXP behavSSF_exact_objective(SEXP theta, SEXP x, SEXP segments,
                                          SEXP k, SEXP context) {
  return objective_impl(theta,x,R_NilValue,R_NilValue,R_NilValue,segments,k,context);
}

extern "C" SEXP behavSSF_simulate(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef methods[] = {
  {"behavSSF_objective", (DL_FUNC) &behavSSF_objective, 7},
  {"behavSSF_exact_objective", (DL_FUNC) &behavSSF_exact_objective, 5},
  {"behavSSF_simulate", (DL_FUNC) &behavSSF_simulate, 7},
  {NULL, NULL, 0}
};
extern "C" void R_init_behavSSF(DllInfo* dll) {
  R_registerRoutines(dll, NULL, methods, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
