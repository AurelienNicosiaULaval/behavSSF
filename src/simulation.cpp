#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Random.h>
#include <R_ext/Utils.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace {
const double pi = 3.14159265358979323846;
struct Node { double x, y; int left, right, axis, point; };
class TargetTree {
  const double* points;
  int count;
  std::vector<int> indices;
  std::vector<Node> nodes;
  int root;
  int build(int lo, int hi, int depth) {
    if (lo >= hi) return -1;
    int mid = lo + (hi-lo)/2, axis = depth%2;
    std::nth_element(indices.begin()+lo, indices.begin()+mid, indices.begin()+hi,
      [&](int a, int b) { return points[a+count*axis] < points[b+count*axis]; });
    int point = indices[mid];
    nodes[mid] = {points[point], points[point+count], -1, -1, axis, point};
    nodes[mid].left = build(lo, mid, depth+1);
    nodes[mid].right = build(mid+1, hi, depth+1);
    return mid;
  }
  void search(int index, double x, double y, double& best_distance, int& best) const {
    if (index < 0) return;
    const Node& node = nodes[index];
    double dx = x-node.x, dy = y-node.y;
    double distance = dx*dx+dy*dy;
    if (distance < best_distance ||
        (distance == best_distance && (best < 0 || node.point < nodes[best].point))) {
      best_distance = distance; best = index;
    }
    double split = node.axis == 0 ? dx : dy;
    int near = split <= 0 ? node.left : node.right;
    int far = split <= 0 ? node.right : node.left;
    search(near, x, y, best_distance, best);
    if (split*split <= best_distance) search(far, x, y, best_distance, best);
  }
public:
  TargetTree(const double* xy, int n): points(xy), count(n), indices(n), nodes(n), root(-1) {
    std::iota(indices.begin(), indices.end(), 0);
    root = build(0, n, 0);
  }
  bool direction(double x, double y, double& bearing) const {
    if (root < 0) return false;
    int best = -1; double distance = R_PosInf;
    search(root, x, y, distance, best);
    if (best < 0 || distance == 0) return false;
    bearing = std::atan2(nodes[best].y-y, nodes[best].x-x);
    return true;
  }
};

int categorical(const double* probabilities, int stride, int k) {
  double u = unif_rand(), sum = 0;
  for (int i = 0; i < k; ++i) {
    sum += probabilities[i*stride];
    if (u <= sum) return i;
  }
  return k-1;
}

// Exact rejection sampler: the uniform proposal has envelope exp(kappa).
// Subtracting kappa keeps the acceptance calculation stable.
double von_mises(double mean, double concentration) {
  for (int attempt = 0; attempt < 1000000; ++attempt) {
    double angle = (2*unif_rand()-1)*pi;
    if (std::log(unif_rand()) <= concentration*(std::cos(angle)-1))
      return std::atan2(std::sin(angle+mean), std::cos(angle+mean));
  }
  throw std::runtime_error("von Mises rejection budget exhausted.");
}
}

extern "C" SEXP behavSSF_simulate(SEXP parameters_, SEXP transition_, SEXP initial_,
                                   SEXP segments_, SEXP targets_, SEXP n_sims_, SEXP fixed_) {
  bool rng_open = false;
  try {
    int k = INTEGER(Rf_getAttrib(parameters_, R_DimSymbol))[0];
    int n_segments = INTEGER(Rf_getAttrib(segments_, R_DimSymbol))[0];
    int n_targets = INTEGER(Rf_getAttrib(targets_, R_DimSymbol))[0];
    int n_sims = Rf_asInteger(n_sims_);
    const double* par = REAL(parameters_);
    const bool hourly_length = INTEGER(Rf_getAttrib(parameters_, R_DimSymbol))[1] == 34;
    const double* transition = REAL(transition_);
    const bool diel = Rf_length(Rf_getAttrib(transition_, R_DimSymbol)) == 3;
    const double* initial = REAL(initial_);
    const double* segments = REAL(segments_);
    const int* fixed = INTEGER(fixed_);
    int total = 0;
    for (int s = 0; s < n_segments; ++s) total += static_cast<int>(segments[s+3*n_segments]);
    bool conditional = Rf_length(fixed_) > 0;
    if (conditional && Rf_length(fixed_) != total)
      throw std::runtime_error("Fixed state path length does not match the segments.");
    TargetTree targets(REAL(targets_), n_targets);
    SEXP out = PROTECT(Rf_allocVector(VECSXP, n_sims));
    for (int rep = 0; rep < n_sims; ++rep) {
      R_CheckUserInterrupt();
      SEXP track = PROTECT(Rf_allocMatrix(REALSXP, total, 5));
      double* result = REAL(track);
      GetRNGstate(); rng_open = true;
      int row = 0;
      for (int s = 0; s < n_segments; ++s) {
        double x = segments[s], y = segments[s+n_segments], heading = segments[s+2*n_segments];
        int length = static_cast<int>(segments[s+3*n_segments]);
        int state = 0;
        double previous_step = 0;
        std::array<double, 26> log_step_history;
        for (int t = 0; t < length; ++t, ++row) {
          const int hour = diel ? (static_cast<int>(segments[s+4*n_segments])+t)%24 : 0;
          if (conditional) state = fixed[row]-1;
          else if (t == 0) state = categorical(initial+hour*k, 1, k);
          else state = categorical(transition+hour*k*k+state, k, k);
          double a = par[state+2*k], b = par[state+3*k];
          double target_bearing = 0;
          if (b != 0 && !targets.direction(x, y, target_bearing)) b = 0;
          double psi = target_bearing-heading;
          double vx = a+b*std::cos(psi), vy = b*std::sin(psi);
          double concentration = std::hypot(vx, vy);
          double angle = von_mises(std::atan2(vy, vx), concentration);
          const double shape = par[state];
          const double rate = par[state+k];
          const double rho = par[state+4*k];
          double scale = 1/rate;
          if (t > 0 && rho != 0) {
            const double log_anchor = std::log(shape)-std::log(rate);
            scale = std::exp(-std::log(rate)+rho*(std::log(previous_step)-log_anchor));
          }
          for (int lag = 22; lag <= 26 && t >= lag; ++lag) {
            const double rho_day = par[state+(5+lag-22)*k];
            if (rho_day != 0) {
              const double log_anchor = std::log(shape)-std::log(rate);
              scale *= std::exp(rho_day*(log_step_history[(t-lag)%26]-log_anchor));
            }
          }
          if (hourly_length) scale *= std::exp(par[state+(10+hour)*k]);
          if (!(scale > 0) || !std::isfinite(scale))
            throw std::runtime_error("Nonpositive or nonfinite gamma scale in simulation.");
          double step = Rf_rgamma(shape, scale);
          if (!(step > 0) || !std::isfinite(step))
            throw std::runtime_error("Nonpositive or nonfinite gamma step generated.");
          heading = std::atan2(std::sin(heading+angle), std::cos(heading+angle));
          x += step*std::cos(heading); y += step*std::sin(heading);
          previous_step = step;
          log_step_history[t%26] = std::log(step);
          result[row] = x; result[row+total] = y; result[row+2*total] = step;
          result[row+3*total] = angle; result[row+4*total] = state+1;
        }
      }
      PutRNGstate(); rng_open = false;
      SET_VECTOR_ELT(out, rep, track);
      UNPROTECT(1);
    }
    UNPROTECT(1);
    return out;
  } catch (const std::exception& e) {
    if (rng_open) PutRNGstate();
    Rf_error("%s", e.what());
  }
  return R_NilValue;
}
