# AGENTS.md

## Purpose

This repository is the optimisation side of a laminate-composite workflow. It is meant to sit above an external FE solver such as CalculiX, not replace one.

The optimiser is responsible for:

1. storing the current design state
2. calling external analysis
3. building local response approximations
4. enforcing laminate feasibility on the optimiser side
5. solving approximate subproblems
6. applying global acceptance and damping logic

Important architectural constraint:

- The FE solver is a response engine and an optional native gradient source
- The FE solver must not be assumed to provide lamination-parameter sensitivities
- lamination-parameter feasibility and as much of the derivative chain as possible stay inside this repo

## Current Cleaned Architecture

The project has been trimmed to the active implementation path.

Active modules:

- `src/OptimisationPipeline/*`
  New orchestration layer and FE-facing interfaces
- `src/GlobalOptimiser/approxFunction.hpp`
  Core response approximation model
- `src/GlobalOptimiser/scminmaxProb.hpp`
  Core separable min-max interior-point solver
- `src/GlobalOptimiser/dampingUtils.hpp`
  Damping utilities used by the outer loop
- `src/Laminate/*`
  Lamination-parameter aggregation and laminate-section feasibility
- `src/Miki/Miki.hpp`
  Miki-cone implementation used by the laminate feasibility layer
- `src/BoundSDP/boundSDP.hpp`
  Scalar lower/upper barrier constraints
- `src/SDPA/SDPA.hpp`
  Fixed-size SDP primitive used by the active cone stack
- `src/SDPA/scalarSDP.hpp`
  Scalar SDP primitive used by `boundSDP`
- `src/SDPA/parameter.hpp`
  Predictor/corrector and step-size controls
- `src/Section/*`
  Dynamic section wrapper used when laminate sections are injected into the core solver

Removed from the active tree:

- `src/InterfaceOptimiser/*`
- `src/GlobalOptimiser/globalOptimiser.*`
- `src/TestProblems/*`
- `src/Example1/*`
- dynamic SDPA refactor files
- generic `fSDPCone` / `LamFSDPCone` branch
- demo-only tests that did not exercise the current architecture

## Source Of Truth

Use these docs in this order when deciding intended behaviour:

1. `docs/gc.pdf`
2. `docs/global.pdf`
3. `docs/sdplam.pdf`
4. `docs/sudo-codes.txt`
5. `docs/Composite-Optimization.drawio.png`

The papers define the intended mathematics. The C++ is an in-progress implementation of that system.

## Code Map

### Orchestration and FE boundary

- `src/OptimisationPipeline/analysis.hpp`
  Analysis request/result types and backend contract
- `src/OptimisationPipeline/approximation.hpp`
  Approximation problem assembly from analysis results
- `src/OptimisationPipeline/subproblem.hpp`
  Generic subproblem interface and simple gradient-penalty fallback solver
- `src/OptimisationPipeline/globalOptimisationDriver.hpp`
  Outer-loop driver with damping, acceptance, finite-difference fallback, and checkpointing
- `src/OptimisationPipeline/jobBackend.hpp`
  Shared file/job backend runner
- `src/OptimisationPipeline/calculixJobBackend.hpp`
  CalculiX-first file/job backend with text-response extraction support
- `src/OptimisationPipeline/defaultAnalysisBackend.hpp`
  Production-facing default backend factory that selects CalculiX by default
- `src/OptimisationPipeline/abaqusJobBackend.hpp`
  Abaqus compatibility wrapper over the shared job backend

### Core subproblem adapters

- `src/OptimisationPipeline/coreMinMaxSubproblem.hpp`
  Adapters from the new pipeline into the core interior-point solver

Current adapters inside that file:

- `CoreMinMax2Var4RespSubproblemSolver`
  Benchmark-shaped bridge for the original 2-variable / 4-response min-max problem
- `CoreLaminateSection1RespSubproblemSolver`
  Direct laminate-aware bridge that injects `laminateSection` objects into the core interior-point solve

### Laminate fallback wrapper

- `src/OptimisationPipeline/laminateAwareSubproblem.hpp`
  Post-solve laminate projection wrapper used as a fallback when the direct core adapter is not the chosen path

### Core numerics

- `src/GlobalOptimiser/approxFunction.hpp`
  Conservative approximation model used by the core solver
- `src/GlobalOptimiser/scminmaxProb.hpp`
  Core predictor/corrector min-max solver
- `src/GlobalOptimiser/dampingUtils.hpp`
  Damping update helpers

Note:

- the directory name `GlobalOptimiser` is historical
- the active contents of that directory are now the core numerical solver pieces, not the old orchestration layer

### Laminate feasibility

- `src/Laminate/varsize.hpp`
- `src/Laminate/cltweight.hpp`
- `src/Laminate/lpfeasible.hpp`
- `src/Laminate/laminateSection.hpp`
- `src/Miki/Miki.hpp`

These files are the optimiser-side laminate feasibility model.

## What The Code Can Do Today

- run the new outer optimisation driver with explicit backend interfaces
- attach FE gradients when available
- generate finite-difference gradients when configured
- solve the benchmark min-max subproblem through the core solver adapter
- solve direct laminate-aware approximate subproblems through the core solver for:
  - one objective
  - optional response constraints
  - multiple laminate section blocks
  - non-zero section offsets
- use a laminate projection wrapper as a fallback path
- run a CalculiX-first job backend with template rendering, command execution, timeout handling, log capture, and result parsing/extraction

## Main Technical Constraints

### 1. The direct laminate core adapter is still narrow

`CoreLaminateSection1RespSubproblemSolver` is the important path, but it is still limited to linearized single-objective laminate problems. It already supports response constraints and multiple section blocks, but it is not yet a fully general production adapter.

### 2. Lamination-parameter sensitivities are still an open production item

The architecture is correct now: solver-native gradients are optional. But the production optimiser-side lamination-parameter derivative path is still unfinished.

### 3. The laminate projection wrapper is now fallback-only

`laminateAwareSubproblem.hpp` is useful, but it should not become the primary laminate integration route if the direct core section-aware path can support the target problem shape.

### 4. The active core solver still lives in historically named modules

The project is cleaner now, but some remaining names are historical:

- `GlobalOptimiser`
- `globopt`
- `scminmaxProb`

Treat these as current core-solver names, not as signs that the old orchestration path still exists.

## Testing

The active test suite is now limited to the current architecture:

- `tests/Miki_test.cpp`
- `tests/calculix_backend_test.cpp`
- `tests/dampingUtils_test.cpp`
- `tests/global_driver_test.cpp`
- `tests/laminate_test.cpp`
- `tests/minmaxprob_test.cpp`
- `tests/section_test.cpp`
- `tests/thickness_test.cpp`

Run with:

- `cmake -S . -B build`
- `cmake --build build -j4`
- `ctest --test-dir build --output-on-failure`

## Recommended Next Step

The next implementation step is:

1. factor laminate section creation/binding out of `CoreLaminateSection1RespSubproblemSolver`
2. make the direct section-aware core solve the default laminate route in orchestration
3. keep `laminateAwareSubproblem.hpp` only as a fallback for unsupported shapes
4. then move to the explicit optimiser-side lamination-parameter derivative path
