# AGENTS.md

## Purpose

This repository is an early-stage research codebase for composite laminate optimisation. It is intended to become an optimisation module that sits above an external finite element analysis workflow, not a standalone structural solver. Abaqus is the clearest target integration, but the architecture should stay general enough to support other FE solvers through the same interface pattern.

Important constraint: Abaqus should be treated primarily as a response engine and only secondarily as an optional native gradient source. The optimisation architecture must not assume that Abaqus can provide sensitivities with respect to lamination parameters. Lamination-parameter mappings, feasibility, and as much of the sensitivity chain as possible need to remain on our side of the interface.

The intended end-to-end system is:

1. Send the current design or laminate-equivalent model to an external FE solver such as Abaqus.
2. Retrieve structural responses and, when available for supported cases, native sensitivities from that solver.
3. Keep the lamination-parameter mapping and laminate feasibility model analytic on our side.
4. Build conservative local approximations of the expensive responses.
5. Solve each approximate subproblem with an interior-point method for separable min-max problems with semidefinite side constraints.
6. Wrap the whole process in a globally convergent outer loop that updates damping until each accepted step is an improvement step.

The project is not finished. The code compiles and the registered tests pass, but several components are still prototypes, partial refactors, or demo programs rather than a unified production solver.

## Source Of Truth

When deciding what the code should do, use these files in this order:

1. `docs/gc.pdf`
   Global outer iteration, improvement-step logic, damping update law.
2. `docs/global.pdf`
   Interior-point method for separable convex min-max problems with PSD side constraints.
3. `docs/sdplam.pdf`
   Lamination-parameter feasibility via SDP and optimisation over the Miki cone.
4. `docs/sudo-codes.txt`
   Fast reference for the algorithms above.
5. `docs/Composite-Optimization.drawio.png`
   High-level architecture picture.

Treat the papers and pseudocode as the intended architecture. Treat the C++ as an incomplete implementation of that architecture.

## Architecture The Project Is Trying To Implement

Conceptually the system is:

`Outer Global Optimiser`
-> owns the optimisation iteration and accepted-design logic
-> builds local approximations around a reference design
-> solves an approximate subproblem
-> evaluates the true responses
-> updates damping
-> accepts only improvement steps

`External Analysis Interface`
-> translates optimisation design variables into FE model inputs
-> launches or drives an external FE solver such as Abaqus
-> extracts responses, constraints, and optional native sensitivities
-> returns that data to the optimisation layer in a solver-agnostic form

`Laminate Parameter Model`
-> keeps the mapping between optimisation variables and lamination parameters analytic
-> enforces lamination-parameter feasibility with Miki-cone based constraints
-> avoids depending on Abaqus for lamination-parameter derivatives

`Approximate Subproblem Solver`
-> min/max reformulation with objective bound `z`
-> interior-point predictor/corrector iterations
-> separable Hessian structure
-> side constraints represented as scalar or matrix SDP cones

`Laminate Feasibility Layer`
-> Miki cone for lamination parameters
-> CLT weights for sublaminate aggregation
-> thickness bounds and laminate section constraints

The intended runtime loop is:

`optimiser / laminate model -> external FE analysis -> approximation builder -> interior-point subproblem solver -> candidate design -> external FE analysis -> acceptance / damping update`

## Code Map

### 1. Global optimisation layer

- `src/GlobalOptimiser/globalOptimiser.hpp`
  Intended outer loop from `gc.pdf`; this should eventually orchestrate FE-solver calls through a clean analysis interface.
- `src/GlobalOptimiser/dampingUtils.hpp`
  Small utility implementation of damping updates.
- `src/GlobalOptimiser/approxFunction.hpp`
  Prototype convex approximation model with sample problems.
- `src/GlobalOptimiser/scminmaxProb.hpp`
  Main legacy implementation of the separable min-max interior-point solver.

### 2. Cone / SDP primitives

- `src/SDPA/parameter.hpp`
  Step-size and predictor/corrector heuristics.
- `src/SDPA/scalarSDP.hpp`
  Scalar barrier constraint building block.
- `src/SDPA/SDPA.hpp`
  Fixed-dimension matrix SDP constraint building block.
- `src/BoundSDP/boundSDP.hpp`
  Lower/upper scalar bound constraints as SDP-compatible side constraints.
- `src/fSDPCone/fSDPCone.hpp`
  Generic linear matrix cone over a basis of matrices.
- `src/Miki/Miki.hpp`
  Miki cone implementations for balanced and unbalanced single-material laminates.

### 3. Laminate-specific aggregation

- `src/Laminate/varsize.hpp`
  Compile-time dimension bookkeeping.
- `src/Laminate/cltweight.hpp`
  CLT aggregation weights across sublaminates.
- `src/Laminate/lpfeasible.hpp`
  Lamination-parameter feasibility solver built by aggregating sublaminate Miki cones.
- `src/Laminate/laminateSection.hpp`
  Adds thickness bounds on top of `lpfeasible`.
- `src/LamFSDPCone/lpfSDPCone.hpp`
  Similar aggregation pattern, but over generic `fSDPCone` building blocks.

### 4. Newer partial abstractions

- `src/SDPA/SDPABase.hpp`
  Dynamic-size abstract interface for SDP solvers.
- `src/SDPA/SDPASolver.hpp`
  Newer generic solver shell, currently still placeholder-level.
- `src/Miki/MikiConeBase.hpp`
  Newer dynamic-size Miki-cone abstraction used by newer tests.

### 5. Legacy wrappers / examples / scaffolding

- `src/InterfaceOptimiser/*`
  Old interface wrappers for response and damping problems; these are the clearest sign that the project is meant to couple to an external analysis code such as Abaqus, but the integration is not complete yet.
- `src/TestProblems/*`
  Small data containers / examples.
- `src/Example1/*`
  C++ inheritance experiments, not part of the optimisation pipeline.

## Current State: What Is Real Vs What Is Draft

### Stable enough to build on

- The repository configures, builds, and runs with:
  - `cmake -S . -B build`
  - `cmake --build build -j4`
  - `ctest --test-dir build --output-on-failure`
- As of `2026-03-23`, all 13 registered tests pass.
- The low-level cone primitives compile and run.
- The laminate feasibility pieces (`Miki`, `lpfeasible`, `laminateSection`) are the most coherent domain-specific part of the code.
- `scminmaxProb.hpp` contains the clearest implementation of the paper-level predictor/corrector logic.
- The overall role of the project is already clear from the docs: this repo is the optimisation engine around an external FE analysis loop, not a replacement for Abaqus or another FE solver.
- The laminate paper already provides the key workaround for the Abaqus sensitivity gap: lamination parameters and their feasibility region are modeled analytically through Miki-cone constructions instead of relying on Abaqus-native lamination-parameter sensitivities.

### Prototype / incomplete / inconsistent parts

- `globopt::GlobalOptimiser::optimize()` is not the current integration point yet.
  The unit test for `globalOptimiser` uses mocks and does not instantiate the full intended path.
- `src/InterfaceOptimiser/responseInterface.hpp` and `src/InterfaceOptimiser/dampingInterface.hpp` do not match the exact interface expected by `src/GlobalOptimiser/globalOptimiser.hpp`.
  The outer-loop integration is therefore not actually complete, especially for a real external FE backend such as Abaqus.
- `src/GlobalOptimiser/scminmaxProb.hpp` still contains multiple comments like `change this once we add lamination`.
  That file is the best current min-max solver implementation, but it is not fully integrated with laminate sections.
- `src/SDPA/SDPASolver.hpp` is not a production solver yet.
  It currently uses placeholder objective data from `GetObjectiveFunctionValues()`.
- `src/Miki/MikiConeBase.hpp` is a simplified dynamic-size refactor, not yet the unified replacement for `src/Miki/Miki.hpp`.
- Several tests under `tests/` are really executable demos with `main()` and prints, not strict correctness tests.
  Passing CTest here mostly means “program exited successfully”.

## Important Technical Seams And Risks

### 1. There are two solver stacks

The codebase currently contains two overlapping architectures:

- Legacy / template-heavy stack:
  `SDPA.hpp`, `scalarSDP.hpp`, `Miki.hpp`, `lpfeasible.hpp`, `laminateSection.hpp`, `scminmaxProb.hpp`
- Newer / dynamic abstraction stack:
  `SDPABase.hpp`, `SDPASolver.hpp`, `MikiConeBase.hpp`

Do not assume these are interchangeable today. They are not. Pick one stack for active development before broad refactoring.

### 2. The papers are more integrated than the code

The intended full pipeline is:

`external FE solve -> response/sensitivity extraction -> local approximation -> min-max interior point solver -> laminate cone side constraints -> outer damping loop`

The repository currently has pieces of that pipeline, but not a single verified executable path that exercises all of them together.

### 3. The FE-coupling boundary should be treated as a first-class design concern

This project should eventually be able to:

- generate solver-ready input from optimisation variables
- run or delegate an Abaqus analysis
- read results back into response objects
- keep the optimisation logic independent from Abaqus-specific file formats as much as possible
- function correctly even when Abaqus provides no native gradients for the active laminate parameterization

That boundary is currently only partially represented in `src/InterfaceOptimiser/*`.

### 4. Do not assume Abaqus can provide lamination-parameter sensitivities

The papers support a different split of responsibilities:

- Abaqus provides high-fidelity responses and, only where supported, native sensitivities for Abaqus-native design variables.
- The laminate model in this repo provides the analytic lamination-parameter feasibility machinery.
- The missing sensitivity path for lamination-parameter optimization should be handled either by:
  - an internal analytic / adjoint-capable laminate-aware solver path, or
  - finite differences in lamination-parameter space as a fallback.

### 5. Test coverage is weaker than the green test suite suggests

Some tests are true assertions, especially the newer GTest-based utility tests. Others only run iterative code and return success if nothing crashes.

Before changing numerics, add tests that check:

- final variable values
- KKT residual reduction
- duality gap reduction
- feasibility preservation
- agreement with known paper examples

### 6. Some “generic” code is not fully generic

- `src/SDPA/SDPA.hpp` uses `Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::Matrix3d>` inside a template.
  That hard-codes a 3x3 assumption in the step-size calculation.
- `src/LamFSDPCone/lpfSDPCone.hpp` exposes `CalulateResiduals` with a typo in the method name.
- Several classes rely on static storage and compile-time dimensions, which makes later runtime-generic integration harder.

## Recommended Development Strategy

### Short-term goal

Do not start with a large refactor. First make the current research pipeline explicit and verifiable.

### Recommended order of work

1. Choose the active implementation path.
   Recommendation: treat the legacy template-heavy path as the current source of executable truth for numerics, and treat the newer `SDPABase` path as exploratory refactor work unless explicitly finishing that refactor.

2. Build one paper-aligned vertical slice.
   A good target is:
   `approxFunction + scminmaxProb + boundSDP`
   first, then extend to
   `laminateSection / lpfeasible`.

3. Define the external analysis contract clearly.
   Before deep integration work, decide what the optimiser expects back from Abaqus or another FE solver:
   mandatory responses, optional native gradients, file-based exchange, process-level API, restart behaviour, and failure handling.

4. Make the lamination-parameter sensitivity path explicit.
   The production architecture should not depend on Abaqus-native lamination-parameter gradients. The preferred paths are:
   - internal analytic / adjoint sensitivity handling around the laminate model
   - finite differences in lamination-parameter space as a fallback

5. Add real regression tests for the papers’ algorithms.
   Convert demo tests into assertion-based tests with tolerances.

6. Integrate the laminate feasibility layer into the min-max solver cleanly.
   This is the biggest missing connection called out in the code comments.

7. Only after the above is stable, refactor toward a single abstraction layer.

## Suggested Near-Term Roadmap

### Milestone 1: Make the current solver verifiable

- Add a small number of deterministic tests around `scminmaxProb`.
- Verify duality gap decreases for the sample problems.
- Verify bound constraints stay feasible.
- Document expected outputs for the existing toy problems.

### Milestone 2: Connect laminate sections to the approximate subproblem solver

- Replace the scalar-side-constraint assumptions in `scminmaxProb.hpp` with a section-based or block-based interface that can host laminate feasibility objects.
- Decide whether a section owns thickness and lamination parameters together, or whether thickness remains a separate side constraint.

### Milestone 3: Rebuild the outer global loop

- Align `GlobalOptimiser`, `responseInterface`, and `dampingInterface`.
- Make the FE-solver interaction explicit and stable, with Abaqus as the first concrete backend and with gradients treated as optional backend capability.
- Make the outer loop accept a real solver object instead of relying on mismatched template conventions.
- Add one end-to-end test that exercises:
  - reference design
  - external analysis call boundary
  - approximation build
  - subproblem solve
  - damping update
  - improvement-step acceptance

### Milestone 4: Build the lamination-parameter sensitivity path outside Abaqus

- Keep lamination-parameter feasibility and LP mappings analytic inside this codebase.
- Decide whether the production sensitivity path is:
  - analytic / adjoint in an internal laminate-aware solver, or
  - finite differences in lamination-parameter space for early production use
- Do not block the architecture on Abaqus-native lamination-parameter sensitivities.

### Milestone 5: Decide whether to keep or replace the newer dynamic abstractions

- If keeping `SDPABase` and `SDPASolver`, finish them and migrate one real cone to that stack.
- Otherwise, remove or isolate them so future work does not split between two architectures.

## Practical Guardrails For Future Agents

- Read the papers first before changing numerical logic.
- Preserve the separation between optimisation logic and FE-backend specifics.
- Preserve the separation between lamination-parameter modeling and FE-backend specifics.
- Prefer extending existing numerical code over “clean rewrite” impulses.
- Keep changes local to one solver stack at a time.
- Add or tighten tests whenever touching:
  - step-size logic
  - residual assembly
  - Hessian assembly
  - damping updates
  - cone feasibility logic
- Preserve Eigen fixed-size types where the existing implementation depends on them for algebra.
- Be careful with sign conventions.
  Many functions implement paper formulas in transformed forms.
- When changing interfaces, verify both the pseudocode intent and the actual call sites.

## Build And Test Commands

Use:

```bash
cmake -S . -B build
cmake --build build -j4
ctest --test-dir build --output-on-failure
```

## First Files To Read Before Any Serious Work

If you are a future agent continuing this project, start here:

1. `docs/sudo-codes.txt`
2. `docs/gc.pdf`
3. `docs/global.pdf`
4. `docs/sdplam.pdf`
5. `src/GlobalOptimiser/scminmaxProb.hpp`
6. `src/Miki/Miki.hpp`
7. `src/Laminate/lpfeasible.hpp`
8. `src/Laminate/laminateSection.hpp`
9. `src/GlobalOptimiser/globalOptimiser.hpp`
10. `tests/minmaxprob_test.cpp`
11. `tests/laminate_test.cpp`

## Bottom Line

This repository already contains much of the mathematical core of the optimisation engine, but not yet a single clean, verified, end-to-end architecture from optimiser to external FE solver. The next useful work is not a cosmetic refactor. It is to make one solver path explicit, tested, paper-faithful, and ready to couple cleanly to Abaqus or another FE backend without assuming that Abaqus will supply lamination-parameter sensitivities.
