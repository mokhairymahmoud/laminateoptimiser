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
- `[x]` Connect the new driver to the core `scminmaxProb` solver as the production subproblem backend.
  - current status: `CoreMinMax2Var4RespSubproblemSolver` covers the benchmark 2-variable / 4-response path
  - current status: `CoreLaminateSection1RespSubproblemSolver` covers direct laminate-aware single-objective linearized solves with response constraints, multiple section blocks, and non-zero offsets
  - current status: `DefaultLaminateSubproblemSolver` now routes supported laminate problems to the direct core solver first and uses laminate projection only as fallback
  - current status: `ConfiguredSubproblemSolver` now provides a production-facing configuration/factory layer for the active benchmark and laminate routing paths

## Phase 4: Laminate-Specific Constraint Integration

- `[x]` Standardize optimiser-side laminate state in the new pipeline request model.
- `[~]` Make laminate sections first-class side constraints in the new subproblem solver path.
  - current status: direct section-aware core solves are in place for the active laminate adapter
  - current status: the direct section-aware route is now the default orchestration choice for supported laminate problems
  - current status: `LaminateAwareSubproblemSolver` remains available as a fallback projection route
  - remaining work: broaden support beyond the current single-objective linearized laminate problem shape
- `[x]` Integrate `laminateSection` ownership cleanly into the orchestration layer.
  - current status: laminate section offsets and bounds flow from `AnalysisRequest` into `ApproximationProblem`
  - current status: the core solver can now carry full `optsection::section` objects with offsets
  - current status: section creation/binding now lives in reusable `laminateSectionBinding` helpers instead of adapter-specific code
  - current status: laminate fallback/projection paths now reuse the same shared laminate section dispatch helpers, so section dispatch is defined in one place
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

- `[x]` Refactor the concrete job-wrapper backend around a solver-neutral file/job runner.
  - current status: common template rendering, command execution, timeout handling, result parsing, and result-file writing now live in a shared job backend layer
  - current status: Abaqus remains available only through a compatibility wrapper
  - current status: a production-facing default backend factory now instantiates CalculiX by default
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
- `[x]` Implement project-specific CalculiX deck generation and real result extraction for the actual FE models used by this project.
- `[x]` Add CalculiX result extraction for the objective and constraint quantities needed by the active approximation pipeline.
  - current status: the repo now contains a concrete CalculiX Composipy benchmark plate deck generator with real composite shell, static displacement, and buckling steps
  - current status: the backend now extracts `buckling_lambda_1` and `tip_u3` from actual CalculiX `.dat` output for that benchmark and assembles benchmark responses around them
  - current status: end-to-end regression coverage now exercises the real benchmark run path when a local `ccx` executable is available
- `[x]` Add a materialized launch path for local CalculiX runs on developer machines:
  - discover or configure the `ccx` executable
  - set job naming and working-directory conventions
  - capture stdout/stderr into run diagnostics
- `[x]` Decide which real CalculiX runs, if any, will provide native sensitivities and which will rely on the optimiser-side/fallback sensitivity path instead.
  - current status: the active Composipy benchmark now declares an explicit response sensitivity policy
  - current status: `mass_per_area` uses backend-native analytic sensitivities
  - current status: `buckling_margin` and `tip_displacement_margin` are explicitly routed through finite-difference fallback for this benchmark
  - current status: the driver now surfaces the configured policy when gradients remain incomplete and fallback is disabled
- `[x]` Keep Abaqus support optional and compatibility-only until there is a concrete need to maintain dual-solver production flows.

## Phase 6: Sensitivities, Fallbacks, and Production Readiness

- `[x]` Support optional FE-supplied sensitivities in `AnalysisResult`.
- `[x]` Add finite-difference fallback sensitivity generation in the driver.
- `[x]` Keep fallback gradients configurable rather than always-on.
- `[x]` Add structured iteration history suitable for diagnostics.
- `[x]` Add basic checkpoint file writing.
- `[x]` Implement the production lamination-parameter sensitivity path outside solver-native gradient support.
  - current status: optimiser-side laminate derivative providers are now plumbed into gradient recovery before finite differences
  - current status: regression coverage now exercises complete-model recovery and mixed backend-plus-laminate gradient assembly on small laminate benchmarks
  - current status: extracted FE scalar quantities can now flow through response-schema assembly and into optimiser-side laminate derivatives via an explicit chain-rule binding layer
  - current status: the active Composipy benchmark now has a production extracted-quantity derivative provider driven by real FE response values, with unit and real-CalculiX regression coverage
- `[x]` Add restart/resume from checkpoint into the active driver flow.
  - current status: checkpoint files now support round-trip parsing back into driver state
  - current status: the active driver can resume from a saved checkpoint design through `optimiseFromCheckpoint(...)`
  - current status: regression coverage verifies resumed runs match an uninterrupted deterministic run
- `[x]` Add richer production logging/output formats.
  - current status: optimisation runs can now be exported as structured `optimisation_summary.json` and `iteration_history.csv` artifacts
  - current status: exported summary data includes final design, objectives, constraints, convergence status, and backend diagnostics
  - current status: regression coverage verifies the structured production artifacts are written for an active driver run
- `[x]` Add failure summaries and diagnostics tailored to real CalculiX runs.
  - current status: failed CalculiX runs now append a backend-specific failure summary with run-directory artifact presence for `job.inp`, `job.dat`, `job.sta`, `job.frd`, and the parsed result file
  - current status: command-failure diagnostics now include sanitized stdout/stderr tails from the captured CalculiX logs
  - current status: invalid-output diagnostics now preserve the extraction error while also reporting which CalculiX artifacts were present or missing
  - current status: regression coverage verifies both command-failure and extraction-failure summaries

## Current Snapshot

- `[x]` The repository has been reduced to the active optimisation architecture.
- `[x]` The build now compiles only active modules and active regression tests.
- `[x]` The active orchestration layer, shared job backend layer, and core solver bridges now build around a CalculiX-first FE wrapper, with Abaqus retained only as a compatibility path.
- `[x]` The core laminate route is in place with shared orchestration-owned section binding helpers across both direct and fallback laminate paths.
- `[x]` Reusable laminate section binding helpers exist in the orchestration layer.
- `[x]` Supported laminate problems now use the direct core laminate route by default, with projection kept as fallback.
- `[x]` The architecture correctly treats solver-native gradients as optional, and the optimiser-side derivative path now supports assembled response terms, partial-gradient merging, an explicit native-vs-fallback benchmark sensitivity policy, and a production extracted-quantity chain for the active Composipy benchmark.

## Next Implementation Step

- `[ ]` Broaden the direct laminate core route beyond the current single-objective linearized laminate problem shape.
