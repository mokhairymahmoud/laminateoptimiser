# End-to-End User Example

This repository is used by wiring three pieces together:

1. an `AnalysisBackend` that runs the external analysis and returns objective/constraint values
2. an `ApproximationBuilder` that linearises the current responses
3. a `SubproblemSolver` that proposes the next design

The smallest runnable example in this repo is:

- [examples/mock_calculix_e2e.cpp](/Users/mohamedkhairy/Documents/projects/comopt/laminateoptimiser/examples/mock_calculix_e2e.cpp)

It uses the active production-style pipeline:

- [src/OptimisationPipeline/calculixJobBackend.hpp](/Users/mohamedkhairy/Documents/projects/comopt/laminateoptimiser/src/OptimisationPipeline/calculixJobBackend.hpp)
- [src/OptimisationPipeline/globalOptimisationDriver.hpp](/Users/mohamedkhairy/Documents/projects/comopt/laminateoptimiser/src/OptimisationPipeline/globalOptimisationDriver.hpp)
- [src/OptimisationPipeline/approximation.hpp](/Users/mohamedkhairy/Documents/projects/comopt/laminateoptimiser/src/OptimisationPipeline/approximation.hpp)
- [src/OptimisationPipeline/subproblem.hpp](/Users/mohamedkhairy/Documents/projects/comopt/laminateoptimiser/src/OptimisationPipeline/subproblem.hpp)

The external solver in the example is a local shell script, not real CalculiX. That keeps it runnable on any machine while preserving the same file-based workflow you would use with CalculiX.

## Build

```bash
cmake -S . -B build
cmake --build build -j4 --target mock_calculix_e2e
```

## Run

```bash
./build/bin/mock_calculix_e2e
```

The example optimises one design variable:

- design variable: plate thickness `t`
- objective: minimise mass, `f(t) = t`
- constraint: enforce displacement limit, `g(t) = 0.18 / t - 0.50 <= 0`

So the optimum is at the active constraint:

- `t ~= 0.36`

## What The Example Does

1. writes a template input file with `{{THICKNESS}}`
2. writes a small mock solver script that reads `job.inp` and produces `job.dat`
3. configures `CalculixJobBackend` to extract:
   - `TOTAL MASS`
   - `TIP DISPLACEMENT`
4. converts those extracted values into:
   - one objective
   - one constraint
5. runs `GlobalOptimisationDriver` with:
   - `LinearApproximationBuilder`
   - `GradientPenaltySubproblemSolver`
   - finite-difference gradient fallback
6. writes:
   - `optimisation_summary.json`
   - `iteration_history.csv`
   - `driver.chk`

## How To Replace The Mock Solver With Real CalculiX

Keep the optimisation code structure the same and change only the backend setup:

1. Replace the generated `plate_template.inp` with your real CalculiX template deck.
2. Point the launch command to `ccx`, or use [src/OptimisationPipeline/defaultAnalysisBackend.hpp](/Users/mohamedkhairy/Documents/projects/comopt/laminateoptimiser/src/OptimisationPipeline/defaultAnalysisBackend.hpp) with `MakeDefaultAnalysisBackend(...)`.
3. Replace the regex extraction rules so they parse the real responses you care about from `.dat`, `.sta`, `.frd`, or any text file your run produces.
4. Keep finite-difference fallback enabled unless your backend or optimiser-side derivative provider supplies the missing gradients.

For a real-CalculiX benchmark already in the repo, see:

- [benchmarks/calculix_e2e_test.cpp](/Users/mohamedkhairy/Documents/projects/comopt/laminateoptimiser/benchmarks/calculix_e2e_test.cpp)
- [benchmarks/calculix_composipy_benchmark_test.cpp](/Users/mohamedkhairy/Documents/projects/comopt/laminateoptimiser/benchmarks/calculix_composipy_benchmark_test.cpp)
