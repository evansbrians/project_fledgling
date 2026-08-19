// rssi_state_hmm.stan ------------------------------------------------------
//
// Two-state hidden Markov model for classifying a radio-tagged bird as
// stationary (state 1) or moving (state 2) from the within-node variance of
// RSSI values, summarised over time windows.
//
// Observation for window n:
//   y[n] = log(pooled within-node RSSI variance)
//
// The model is fit to one tag at a time (see detect_movement_states() in
// rssi_state_hmm.R, which maps this model across all tags). State means are
// ordered by construction (state 2 is always the higher-variance state), so
// labels cannot switch between chains or tags.
//
// The stationary state is defined from the full population of tags: the
// 22:00-03:00 reference period enters twice:
//   (1) as the prior location for the stationary emission mean
//       (mu_ref_mean / mu_ref_sd, computed in R from the night windows of
//       ALL tags, not just the tag being fit), and
//   (2) as a soft label on this tag's night windows: their emission
//       densities are tilted toward the stationary state by ref_leak, the
//       prior probability that a bird is moving during the reference
//       period. Setting ref_leak = 0 would hard-label night windows; a
//       small non-zero value (~0.05) lets genuinely active nights override
//       it.

functions {
  
  // K x N matrix of log observation densities, with the night-window tilt
  // folded in.
  
  matrix rssi_log_omega(vector y,
                        array[] int is_ref,
                        vector mu,
                        vector sigma,
                        real nu,
                        real ref_leak) {
    int n = num_elements(y);
    matrix[2, n] log_omega;

    for (i in 1:n) {
      log_omega[1, i] = student_t_lpdf(y[i] | nu, mu[1], sigma[1]);
      log_omega[2, i] = student_t_lpdf(y[i] | nu, mu[2], sigma[2]);

      if (is_ref[i] == 1) {
        log_omega[1, i] += log1m(ref_leak);
        log_omega[2, i] += log(ref_leak);
      }
    }

    return log_omega;
  }
}

data {
  int<lower = 1> n_win;                             // windows, one tag
  int<lower = 1> n_seg;                             // gap-free runs of windows

  vector[n_win] y;                                  // log pooled RSSI variance
  array[n_win] int<lower = 0, upper = 1> is_ref;    // night reference window

  array[n_seg] int<lower = 1, upper = n_win> seg_start;
  array[n_seg] int<lower = 1> seg_len;

  // Population reference-period prior for the stationary emission mean:
  
  real mu_ref_mean;
  real<lower = 0> mu_ref_sd;

  // Prior on log(separation between state means):
  
  real delta_mean;
  real<lower = 0> delta_sd;

  array[2] vector<lower = 0>[2] alpha_gamma;        // Dirichlet, transitions
  real<lower = 0, upper = 0.5> ref_leak;
  real<lower = 1> nu;                               // Student-t d.f.

  int<lower = 0, upper = 1> save_states;            // 1 = return p_moving
}

parameters {
  real mu_stat;                                     // stationary emission mean
  real log_delta;                                   // log separation of means
  vector<lower = 0>[2] sigma;                       // emission scale by state

  array[2] simplex[2] gamma_rows;                   // transition matrix rows
  simplex[2] rho;                                   // initial state
}

transformed parameters {
  vector[2] mu;

  mu[1] = mu_stat;
  mu[2] = mu_stat + exp(log_delta);
}

model {
  matrix[2, 2] gamma_mat;
  matrix[2, n_win] log_omega;

  mu_stat ~ normal(mu_ref_mean, mu_ref_sd);
  log_delta ~ normal(delta_mean, delta_sd);
  sigma ~ student_t(3, 0, 1);

  for (k in 1:2) {
    gamma_rows[k] ~ dirichlet(alpha_gamma[k]);
    gamma_mat[k] = gamma_rows[k]';
  }

  log_omega = rssi_log_omega(y, is_ref, mu, sigma, nu, ref_leak);

  // One forward pass per gap-free run of windows, so that a break in
  // detections is never modelled as a state transition.
  
  for (s in 1:n_seg) {
    target += hmm_marginal(
      block(log_omega, 1, seg_start[s], 2, seg_len[s]),
      gamma_mat,
      rho
    );
  }
}

generated quantities {
  real p_stay_stationary = gamma_rows[1, 1];
  real p_stay_moving = gamma_rows[2, 2];
  vector[save_states ? n_win : 0] p_moving;

  if (save_states) {
    matrix[2, 2] gamma_mat;
    matrix[2, n_win] log_omega;

    for (k in 1:2) {
      gamma_mat[k] = gamma_rows[k]';
    }

    log_omega = rssi_log_omega(y, is_ref, mu, sigma, nu, ref_leak);

    for (s in 1:n_seg) {
      matrix[2, seg_len[s]] state_prob =
        hmm_hidden_state_prob(
          block(log_omega, 1, seg_start[s], 2, seg_len[s]),
          gamma_mat,
          rho
        );

      for (i in 1:seg_len[s]) {
        p_moving[seg_start[s] + i - 1] = state_prob[2, i];
      }
    }
  }
}
