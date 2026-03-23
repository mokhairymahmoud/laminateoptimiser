# Laminate Optimiser Implementation Plan

This file tracks the implementation roadmap for the cleaned project structure.

Important constraint:

- The external FE solver is a response engine and optional native gradient source
- The FE solver is not assumed to provide lamination-parameter sensitivities

Status legend:

- `[x]` implemented
- `[~]` partially implemented
- `[ ]` not implemented yet

## Phase 0: Repository Cleanup

- `[x]` Remove abandoned interface/demo branches from the active project:
  - old `InterfaceOptimiser` wrappers
  - old `GlobalOptimiser/globalOptimiser.*` orchestration shell
  - `TestProblems`
  - `Example1`
  - dynamic SDPA refactor branch
  - generic `fSDPCone` / `LamFSDPCone` branch
- `[x]` Switch CMake from broad globs to explicit active source and test lists.
- `[x]` Rename the active pipeline-to-core adapters from `legacy*` to `core*`.
- `[x]` Replace demo-only section testing with a focused regression on the active dynamic section dispatch path.

## Phase 1: Numerical Core Hardening

- `[x]` Treat `scminmaxProb + boundSDP + Miki + lpfeasible + laminateSection` as the active numerical path.
- `[x]` Remove non-active solver branches from the cleaned repository tree.
- `[x]` Convert the kept core benchmarks into assertion-based numerical regressions for:
  - `minmaxprob`
  - `Miki`
  - `thickness`
  - `laminate feasibility`
  - dynamic laminate section dispatch
- `[x]` Define expected regression outputs for benchmark problems.
- `[x]` Verify the cleaned suite with `ctest --test-dir build --output-on-failure`.
- `[ ]` Audit and fix remaining blocking numerical inconsistencies in the core solver path beyond the currently covered regressions.

## Phase 2: Explicit Optimisation Pipeline Interfaces

- `[x]` Add explicit pipeline types:
  - `AnalysisRequest`
  - `AnalysisResult`
  - `AnalysisBackend`
  - `ApproximationBuilder`
  - `SubproblemSolver`
- `[x]` Add the orchestration module under `src/OptimisationPipeline/`.
- `[x]` Support optional gradients in `AnalysisResult`.
- `[x]` Keep backend gradients optional rather than mandatory.
- `[x]` Remove the abandoned old orchestration/interface path from the repository.

## Phase 3: End-to-End Optimisation Driver

- `[x]` Add `GlobalOptimisationDriver`.
- `[x]` Implement outer-loop acceptance/rejection behavior with damping updates.
- `[x]` Keep FE-backend interaction separate from damping and acceptance logic.
- `[x]` Add deterministic stopping conditions for:
  - outer-iteration limit
  - sub-iteration limit
  - stagnation
  - backend failure
- `[x]` Keep serializable run state/history and add checkpoint writing support.
- `[x]` Add integration tests with fake backends and scripted subproblem solvers.
- `[~]` Connect the new driver to the core `scminmaxProb` solver as the production subproblem backend.
  - current status: `CoreMinMax2Var4RespSubproblemSolver` covers the benchmark 2-variable / 4-response path
  - current status: `CoreLaminateSection1RespSubproblemSolver` covers direct laminate-aware single-objective linearized solves with response constraints, multiple section blocks, and non-zero offsets
  - current status: `DefaultLaminateSubproblemSolver` now routes supported laminate problems to the direct core solver first and uses laminate projection only as fallback
  - remaining work: consolidate the specialized adapters behind a cleaner production-facing configuration/factory layer

## Phase 4: Laminate-Specific Constraint Integration

- `[x]` Standardize optimiser-side laminate state in the new pipeline request model.
- `[~]` Make laminate sections first-class side constraints in the new subproblem solver path.
  - current status: direct section-aware core solves are in place for the active laminate adapter
  - current status: the direct section-aware route is now the default orchestration choice for supported laminate problems
  - current status: `LaminateAwareSubproblemSolver` remains available as a fallback projection route
  - remaining work: broaden support beyond the current single-objective linearized laminate problem shape
- `[~]` Integrate `laminateSection` ownership cleanly into the orchestration layer.
  - current status: laminate section offsets and bounds flow from `AnalysisRequest` into `ApproximationProblem`
  - current status: the core solver can now carry full `optsection::section` objects with offsets
  - current status: section creation/binding now lives in reusable `laminateSectionBinding` helpers instead of adapter-specific code
  - remaining work: reuse the same helpers in every laminate fallback/projection path so section dispatch is defined in one place
- `[~]` Assemble approximate subproblems that simultaneously include:
  - response constraints
  - laminate feasibility constraints
  - bound constraints
  - current status: the direct core laminate adapter already does this for the cleaned single-objective linearized laminate path
  - remaining work: broaden that path beyond the current objective/approximation shape
- `[~]` Make the lamination-parameter mapping and derivative path explicit on our side of the architecture rather than relying on solver-native gradients.
  - current status: optimiser-owned lamination-parameter derivative interfaces now exist in the pipeline
  - current status: assembled laminate response terms can now mix dense design contributions with laminate-section contributions
  - current status: the driver now merges masked backend gradients, optimiser-side laminate gradients, and finite-difference fallback in one recovery path
  - remaining work: connect the assembled response terms to the real FE response quantities extracted for active runs

## Phase 5: CalculiX-First FE Job-Wrapper Backend

- `[~]` Refactor the concrete job-wrapper backend around a solver-neutral file/job runner.
  - current status: common template rendering, command execution, timeout handling, and result parsing now live in a shared job backend layer
  - current status: Abaqus remains available through a compatibility wrapper while the integration is being migrated
  - current status: a CalculiX-specific backend wrapper now exists and is the intended active FE entry point
  - remaining work: move the orchestration and project configuration to instantiate CalculiX by default
- `[x]` Support template input rendering from optimisation variables.
- `[x]` Support isolated run directories.
- `[x]` Support external command execution.
- `[x]` Support timeout and failure reporting.
- `[x]` Support parsing result files into `AnalysisResult`.
- `[x]` Treat solver-native gradients as optional backend output, not as a required capability.
- `[x]` Add backend tests for:
  - parameter injection
  - result parsing
  - command failure
  - timeout handling
- `[ ]` Implement project-specific CalculiX deck generation and real result extraction for the actual FE models used by this project.
- `[ ]` Add CalculiX result extraction for the objective and constraint quantities needed by the active approximation pipeline.
- `[ ]` Add a materialized launch path for local CalculiX runs on developer machines:
  - discover or configure the `ccx` executable
  - set job naming and working-directory conventions
  - capture stdout/stderr into run diagnostics
- `[ ]` Decide which real CalculiX runs, if any, will provide native sensitivities and which will rely on the optimiser-side/fallback sensitivity path instead.
- `[ ]` Keep Abaqus support optional and compatibility-only until there is a concrete need to maintain dual-solver production flows.

## Phase 6: Sensitivities, Fallbacks, and Production Readiness

- `[x]` Support optional FE-supplied sensitivities in `AnalysisResult`.
- `[x]` Add finite-difference fallback sensitivity generation in the driver.
- `[x]` Keep fallback gradients configurable rather than always-on.
- `[x]` Add structured iteration history suitable for diagnostics.
- `[x]` Add basic checkpoint file writing.
- `[~]` Implement the production lamination-parameter sensitivity path outside solver-native gradient support.
  - current status: optimiser-side laminate derivative providers are now plumbed into gradient recovery before finite differences
  - current status: regression coverage now exercises complete-model recovery and mixed backend-plus-laminate gradient assembly on small laminate benchmarks
  - remaining work: replace the in-test response assembly with the production laminate derivative chain driven by real FE response extraction
- `[ ]` Add restart/resume from checkpoint into the active driver flow.
- `[ ]` Add richer production logging/output formats.
- `[ ]` Add failure summaries and diagnostics tailored to real CalculiX runs.

## Current Snapshot

- `[x]` The repository has been reduced to the active optimisation architecture.
- `[x]` The build now compiles only active modules and active regression tests.
- `[~]` The active orchestration layer, shared job backend layer, and core solver bridges build together while the FE wrapper is being migrated from Abaqus-first to CalculiX-first.
- `[~]` The core laminate route is in place with shared orchestration-owned section binding helpers, but fallback/projection paths do not all reuse the same helper layer yet.
- `[x]` Reusable laminate section binding helpers exist in the orchestration layer.
- `[x]` Supported laminate problems now use the direct core laminate route by default, with projection kept as fallback.
- `[~]` The architecture correctly treats solver-native gradients as optional, and the optimiser-side lamination-parameter derivative path now supports assembled response terms plus partial-gradient merging, but the production derivative chain is still unfinished.

## Next Implementation Step

- `[ ]` Make the CalculiX backend the default FE-facing path and bind the assembled optimiser-side laminate response terms to the real result extraction path.

Concretely:

1. define the CalculiX run contract:
   - input deck templating conventions
   - result-file locations
   - objective/constraint extraction format
2. switch orchestration/configuration to construct CalculiX-backed analysis runs by default
3. build the production response-term assembly from CalculiX-extracted quantities instead of constructing benchmark terms in tests
4. keep masked backend/native gradients, optimiser-side laminate terms, and finite-difference fallback interoperable for unsupported pieces
5. add regression coverage around the real CalculiX-facing extraction and assembly path once that contract is defined
