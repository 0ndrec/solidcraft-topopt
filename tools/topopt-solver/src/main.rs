//! TopOpt Solver: topology optimization on all CPU cores. The command line is
//! [`topopt_solver::cli`]; see docs/topopt-solver.md.

fn main() {
    topopt_solver::cli::main_with(topopt_solver::cli::App::cpu(), std::env::args().skip(1).collect())
}
