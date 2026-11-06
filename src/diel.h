#ifndef BEHAVSSF_DIEL_H
#define BEHAVSSF_DIEL_H

// Transition from the preceding state to the state at the specified hour.
// There are two logit regressions, one for each direction of switching.
class DielTransitions {
  int q;
  const double* basis;
  std::vector<double> probability, log_probability, initial, initial_gradient;

public:
  DielTransitions(const double* basis_, int q_, const double* coefficients):
    q(q_), basis(basis_), probability(24*4), log_probability(24*4),
    initial(24*2), initial_gradient(24*2*q) {
    double u = 0, v = 1;
    std::vector<double> du(2*q, 0), dv(2*q, 0);
    for (int h = 0; h < 24; ++h) {
      for (int i = 0; i < 2; ++i) {
        double eta = 0;
        for (int j = 0; j < q; ++j) eta += basis[h+24*j]*coefficients[j+q*i];
        if (!std::isfinite(eta)) throw std::runtime_error("Nonfinite diel transition logit.");
        const double z = std::log1p(std::exp(-std::abs(eta)));
        const double stay = -std::max(0.0, eta)-z;
        const double move = -std::max(0.0, -eta)-z;
        log_probability[4*h+3*i] = stay;
        log_probability[4*h+1+i] = move;
        probability[4*h+3*i] = std::exp(stay);
        probability[4*h+1+i] = std::exp(move);
      }
      const double* p = at(h);
      const double a = p[1], b = p[2], c = 1-a-b;
      // Compose q_next = a + c*q_prev and differentiate the composition.
      for (int j = 0; j < 2*q; ++j) {
        const double x = basis[h+24*(j%q)];
        const double da = j < q ? x*a*p[0] : 0;
        const double db = j >= q ? x*b*p[3] : 0;
        du[j] = da - (da+db)*u + c*du[j];
        dv[j] = -(da+db)*v + c*dv[j];
      }
      u = a+c*u;
      v = c*v;
    }
    // q at the end of the day must equal q before the first hour.
    double long_state = u/(1-v);
    std::vector<double> derivative(2*q);
    for (int j = 0; j < 2*q; ++j) derivative[j] = (du[j]+long_state*dv[j])/(1-v);
    for (int h = 0; h < 24; ++h) {
      const double* p = at(h);
      const double a = p[1], b = p[2], c = 1-a-b;
      for (int j = 0; j < 2*q; ++j) {
        const double x = basis[h+24*(j%q)];
        const double da = j < q ? x*a*p[0] : 0;
        const double db = j >= q ? x*b*p[3] : 0;
        derivative[j] = da*(1-long_state) - db*long_state + c*derivative[j];
        initial_gradient[2*q*h+j] = derivative[j];
      }
      long_state = a+c*long_state;
      if (!(long_state > 0 && long_state < 1))
        throw std::runtime_error("Invalid periodic initial probabilities.");
      initial[2*h] = 1-long_state;
      initial[2*h+1] = long_state;
    }
  }

  const double* at(int h) const { return probability.data()+4*h; }
  const double* log_at(int h) const { return log_probability.data()+4*h; }
  const double* start_at(int h) const { return initial.data()+2*h; }

  void add_start_score(int h, const std::vector<double>& posterior, double* score) const {
    const double factor = posterior[1]/initial[2*h+1] - posterior[0]/initial[2*h];
    for (int j = 0; j < 2*q; ++j) score[j] += factor*initial_gradient[2*q*h+j];
  }

  void add_transition_score(int h, int i, double residual, double* score) const {
    for (int j = 0; j < q; ++j) score[j+q*i] += residual*basis[h+24*j];
  }
};

#endif
