#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

// TODO: Set the timestep length and duration [DONE]

// ======================= MY IMPLEMENTATION | START ======================= //

// T = 10 * 0.1 = 1 seconds (can tune up later)
size_t N_TIME_STEPS = 10; // unsigned integer
double dt = 0.1;
// For setting numeric limits
// https://www.npr.org/sections/krulwich/2012/09/17/161096233/which-is-greater-the-number-of-sand-grains-on-earth-or-stars-in-the-sky
const double SAND_PARTICLES_ON_EARTH = 7.5e18;
// Steering angle limit, pre-calculate 25 degrees to radians
const double STEERING_MAX_LIMIT = 0.436332;
// Throttle limit
const double THROTTLE_MAX_LIMIT = 1.0;
// State variable and actuation variables
const size_t N_STATE_VARS = 6;
const size_t N_ACTUATION_VARS = 2;

// ======================== MY IMPLEMENTATION | END ======================== //

// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
// presented in the classroom matched the previous radius.
//
// This is the length from front to CoG that has a similar radius.
const double Lf = 2.67;

// ======================= MY IMPLEMENTATION | START ======================= //

// Reference velocity is used to penalize for driving too slow or too fast
const double ref_v = 35; // can tune later

// Initialize indices for convenient access into `fg` vector of cost constraints
size_t x_start_idx = 0;
size_t y_start_idx = x_start_idx + N_TIME_STEPS;
size_t psi_start_idx = y_start_idx + N_TIME_STEPS;
size_t v_start_idx = psi_start_idx + N_TIME_STEPS;
size_t cte_start_idx = v_start_idx + N_TIME_STEPS;
size_t psi_err_start_idx = cte_start_idx + N_TIME_STEPS;
size_t delta_start_idx = psi_err_start_idx + N_TIME_STEPS;
size_t a_start_idx = delta_start_idx + N_TIME_STEPS - 1;

// ======================== MY IMPLEMENTATION | END ======================== //

class FG_eval {
 public:
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;
  FG_eval(Eigen::VectorXd coeffs) { this->coeffs = coeffs; }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars) {
    // TODO: implement MPC [DONE]
    // `fg` a vector of the cost constraints, `vars` is a vector of variable values (state & actuators)
    // NOTE: You'll probably go back and forth between this function and
    // the Solver function below.

// ======================= MY IMPLEMENTATION | START ======================= //

    // ----- Calculate cost | start ----- //

    // Initilize cost to be 0 and then update
    fg[0] = 0;

    // The part of the cost based on the reference state
    for (int t = 0; t < N_TIME_STEPS; t++) {
      // Penalize for cross-track error
      fg[0] += CppAD::pow(vars[cte_start_idx + t], 2);
      // Penalize for orientation angle
      fg[0] += CppAD::pow(vars[psi_err_start_idx + t], 2);
      // Penalize for stopping or driving too fast
      fg[0] += CppAD::pow(vars[v_start_idx + t] - ref_v, 2);
    }

    // Minimize the use of steering and throttle / brake
    for (int t = 0; t < N_TIME_STEPS - 1; t++) {
      fg[0] += CppAD::pow(vars[delta_start_idx + t], 2);
      fg[0] += CppAD::pow(vars[a_start_idx + t], 2);
    }

    // Minimize fast changing of steering and throttle / brake
    for (int t = 0; t < N_TIME_STEPS - 2; t++) {
      fg[0] += CppAD::pow(vars[delta_start_idx + t + 1] - vars[delta_start_idx + t], 2);
      fg[0] += CppAD::pow(vars[a_start_idx + t + 1] - vars[a_start_idx + t], 2);
    }
    // ------ Calculate cost | end ------ //


    // ----- Calculate constraints | start ----- //

    // Variable at position `t` in `vars`
    // is stored at `t+1` in `fg` because fg[0] is cost
    fg[1 + x_start_idx] = vars[x_start_idx];
    fg[1 + y_start_idx] = vars[y_start_idx];
    fg[1 + psi_start_idx] = vars[psi_start_idx];
    fg[1 + v_start_idx] = vars[v_start_idx];
    fg[1 + cte_start_idx] = vars[cte_start_idx];
    fg[1 + psi_err_start_idx] = vars[psi_err_start_idx];

    // The rest of the constraints
    for (int t = 1; t < N_TIME_STEPS; t++) {
      // The state at time t+1
      // Retreive variables from vectors for convenience
      AD<double> x1 = vars[x_start_idx + t];
      AD<double> y1 = vars[y_start_idx + t];
      AD<double> psi1 = vars[psi_start_idx + t];
      AD<double> v1 = vars[v_start_idx + t];
      AD<double> cte1 = vars[cte_start_idx + t];
      AD<double> psi_err1 = vars[psi_err_start_idx + t];

      // The state at time t
      // Retreive variables from vectors for convenience
      AD<double> x0 = vars[x_start_idx + t - 1];
      AD<double> y0 = vars[y_start_idx + t - 1];
      AD<double> psi0 = vars[psi_start_idx + t - 1];
      AD<double> v0 = vars[v_start_idx + t - 1];
      AD<double> cte0 = vars[cte_start_idx + t - 1];
      AD<double> psi_err0 = vars[psi_err_start_idx + t - 1];

      // Only consider the actuation at time t
      AD<double> delta0 = vars[delta_start_idx + t - 1];
      AD<double> a0 = vars[a_start_idx + t - 1];

      // Target path polynomial, calculate desired orientation
      // 1. Pre-calculate for convenience
      AD<double> x0_square = pow(x0, 2);
      AD<double> x0_cube = pow(x0, 3);
      // 2. Use 3rd degree polynomial
      AD<double> f0 = coeffs[0] + coeffs[1] * x0 + coeffs[2] * x0_square + coeffs[3] * x0_cube;
      // 3. Take atan of derivative of 3rd degree polynomial
      AD<double> f0_prime = coeffs[1] + 2 * coeffs[2] * x0 + 3 * coeffs[3] * x0_square;
      AD<double> psi_des0 = CppAD::atan(f0_prime);

      // Calculate next state according to kinematic model: [x, y, psi, v]
      fg[1 + x_start_idx + t] = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);
      fg[1 + y_start_idx + t] = y1 - (y0 + v0 * CppAD::sin(psi0) * dt);
      fg[1 + psi_start_idx + t] = psi1 - (psi0 + v0 * delta0 / Lf * dt);
      fg[1 + v_start_idx + t] = v1 - (v0 + a0 * dt);

      // Calculate next cross-track error and orientation error: [cte, psi_err]
      fg[1 + cte_start_idx + t] = cte1 - ((f0 - y0) + (v0 * CppAD::sin(psi_err0) * dt));
      fg[1 + psi_err_start_idx + t] = psi_err1 - ((psi0 - psi_des0) + v0 * delta0 / Lf * dt);
    }

    // ---- Calculate constraints | end ---- //

// ======================== MY IMPLEMENTATION | END ======================== //

  }
};

//
// MPC class definition implementation.
//
MPC::MPC() {}
MPC::~MPC() {}

vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {
  bool ok = true;
  size_t i;
  typedef CPPAD_TESTVECTOR(double) Dvector;

  // TODO: Set the number of model variables (includes both states and inputs) [DONE]

// ======================= MY IMPLEMENTATION | START ======================= //

  // Six state variables, each will have N_TIME_STEPS time steps
  // Two actuation variables, each will have N_TIME_STEPS - 1 steps

  // Total number of vars
  const size_t n_vars = N_STATE_VARS * N_TIME_STEPS + N_ACTUATION_VARS * (N_TIME_STEPS - 1);

  // For convenience, extract state variables
  double x = state[0];
  double y = state[1];
  double psi = state[2];
  double v = state[3];
  double cte = state[4];
  double psi_err = state[5];

// ======================== MY IMPLEMENTATION | END ======================== //

  // TODO: Set the number of constraints [DONE]

// ======================= MY IMPLEMENTATION | START ======================= //

  // One constraint for each state variable
  size_t n_constraints = N_STATE_VARS * N_TIME_STEPS;

// ======================== MY IMPLEMENTATION | END ======================== //

  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  Dvector vars(n_vars);
  for (int i = 0; i < n_vars; i++) {
    vars[i] = 0;
  }

  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);
  // TODO: Set lower and upper limits for variables [DONE]

// ======================= MY IMPLEMENTATION | START ======================= //

  // Need to set everything except for steering and throttle to large limits
  for (int j = 0; j < delta_start_idx; ++j) {
    // How about estimated number of sand particles on earth?
    vars_lowerbound[j] = -SAND_PARTICLES_ON_EARTH;
    vars_upperbound[j] = SAND_PARTICLES_ON_EARTH;
  }

  // Steering limits in radians (pre-calculated)
  for (int k = delta_start_idx; k < a_start_idx; ++k) {
    vars_lowerbound[k] = -STEERING_MAX_LIMIT;
    vars_upperbound[k] = STEERING_MAX_LIMIT;
  }

  // Throttle limits
  for (int l = a_start_idx; l < n_vars; ++l) {
    vars_lowerbound[l] = -THROTTLE_MAX_LIMIT;
    vars_upperbound[l] = THROTTLE_MAX_LIMIT;
  }

// ======================== MY IMPLEMENTATION | END ======================== //

  // Lower and upper limits for the constraints
  // Should be 0 besides initial state.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (int i = 0; i < n_constraints; i++) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }

// ======================= MY IMPLEMENTATION | START ======================= //

  // Also need to set constraints lower and upper bounds to current values of state
  constraints_lowerbound[x_start_idx] = x;
  constraints_lowerbound[y_start_idx] = y;
  constraints_lowerbound[psi_start_idx] = psi;
  constraints_lowerbound[v_start_idx] = v;
  constraints_lowerbound[cte_start_idx] = cte;
  constraints_lowerbound[psi_err_start_idx] = psi_err;

  constraints_upperbound[x_start_idx] = x;
  constraints_upperbound[y_start_idx] = y;
  constraints_upperbound[psi_start_idx] = psi;
  constraints_upperbound[v_start_idx] = v;
  constraints_upperbound[cte_start_idx] = cte;
  constraints_upperbound[psi_err_start_idx] = psi_err;

// ======================== MY IMPLEMENTATION | END ======================== //

  // object that computes objective and constraints
  FG_eval fg_eval(coeffs);

  //
  // NOTE: You don't have to worry about these options
  //
  // options for IPOPT solver
  std::string options;
  // Uncomment this if you'd like more print information
  options += "Integer print_level  0\n";
  // NOTE: Setting sparse to true allows the solver to take advantage
  // of sparse routines, this makes the computation MUCH FASTER. If you
  // can uncomment 1 of these and see if it makes a difference or not but
  // if you uncomment both the computation time should go up in orders of
  // magnitude.
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit.
  options += "Numeric max_cpu_time          0.5\n";

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  // Cost
  auto cost = solution.obj_value;
  std::cout << "Cost " << cost << std::endl;

  // TODO: Return the first actuator values. The variables can be accessed with [DONE]
  // `solution.x[i]`.
  //
  // {...} is shorthand for creating a vector, so auto x1 = {1.0,2.0}
  // creates a 2 element double vector.

// ======================= MY IMPLEMENTATION | START ======================= //

  vector<double> result = {solution.x[delta_start_idx], solution.x[a_start_idx]};
  for (int m = 0; m < N_TIME_STEPS; ++m) {
    result.push_back(solution.x[2 + m]);
  }
  for (int n = 0; n < N_TIME_STEPS; ++n) {
    result.push_back(solution.x[2 + N_TIME_STEPS + n]);
  }
  return result;

// ======================== MY IMPLEMENTATION | END ======================== //
}
