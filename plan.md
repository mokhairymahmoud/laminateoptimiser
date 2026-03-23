# Laminate Optimiser Implementation Plan

This file tracks the implementation roadmap for turning the project into a laminate optimisation module that can couple to external FE solvers such as Abaqus.

Important constraint: Abaqus is treated as a response engine and optional native gradient source, not as the assumed provider of lamination-parameter sensitivities. The production architecture must support lamination-parameter optimization even when Abaqus cannot return those gradients directly.

Status legend:

- `[x]` implemented
- `[~]` partially implemented
- `[ ]` not implemented yet

## Phase 1: Numerical Core Hardening

- `[x]` Treat `scminmaxProb + boundSDP + Miki + lpfeasible + laminateSection` as the active numerical path.
- `[x]` Keep the newer dynamic solver stack non-authoritative for now.
- `[x]` Convert demo-style tests into assertion-based numerical regressions for:
  - `minmaxprob`
  - `Miki`
  - `thickness`
  - `laminate feasibility`
  - `fSDPCone`
- `[x]` Define expected regression outputs for benchmark problems.
- `[x]` Verify the suite with `ctest --test-dir build --output-on-failure`.
- `[ ]` Audit and fix remaining blocking numerical inconsistencies in the legacy solver path beyond the currently covered regressions.

## Phase 2: Explicit Optimisation Pipeline Interfaces

- `[x]` Add explicit pipeline types:
  - `AnalysisRequest`
  - `AnalysisResult`
  - `AnalysisBackend`
  - `ApproximationBuilder`
  - `SubproblemSolver`
- `[x]` Add a new orchestration module under `src/OptimisationPipeline/`.
- `[x]` Support optional gradients in `AnalysisResult`.
- `[x]` Keep backend gradients optional rather than mandatory.
- `[x]` Preserve the FE-coupling boundary as a dedicated interface instead of extending the legacy wrappers ad hoc.
- `[ ]` Replace or fully wrap the old `GlobalOptimiser` / `InterfaceOptimiser` path behind the new interfaces.

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
- `[ ]` Connect the new driver to the legacy `scminmaxProb` solver as the production subproblem backend.

## Phase 4: Laminate-Specific Constraint Integration

- `[x]` Standardize optimiser-side laminate state in the new pipeline request model.
- `[ ]` Make the lamination-parameter mapping and derivative path explicit on our side of the architecture rather than relying on Abaqus-native gradients.
- `[ ]` Make laminate sections first-class side constraints in the new subproblem solver path.
- `[ ]` Integrate `laminateSection` ownership cleanly into the new orchestration layer.
- `[ ]` Assemble approximate subproblems that simultaneously include:
  - response constraints
  - laminate feasibility constraints
  - bound constraints
- `[ ]` Validate the new pipeline using laminate-aware subproblem solves rather than only standalone laminate regressions.

## Phase 5: Abaqus Job-Wrapper Backend

- `[x]` Add a concrete Abaqus-style job-wrapper backend.
- `[x]` Support template input rendering from optimisation variables.
- `[x]` Support isolated run directories.
- `[x]` Support external command execution.
- `[x]` Support timeout and failure reporting.
- `[x]` Support parsing result files into `AnalysisResult`.
- `[x]` Treat Abaqus-native gradients as optional backend output, not as a required capability.
- `[x]` Add backend tests for:
  - parameter injection
  - result parsing
  - command failure
  - timeout handling
- `[ ]` Implement project-specific Abaqus deck generation and real result extraction for the actual FE models used by this project.
- `[ ]` Decide which real Abaqus runs, if any, will provide native sensitivities and which will rely on the internal / fallback sensitivity path instead.

## Phase 6: Sensitivities, Fallbacks, and Production Readiness

- `[x]` Support optional FE-supplied sensitivities in `AnalysisResult`.
- `[x]` Add finite-difference fallback sensitivity generation in the driver.
- `[x]` Keep fallback gradients configurable rather than always-on.
- `[x]` Add structured iteration history suitable for diagnostics.
- `[x]` Add basic checkpoint file writing.
- `[ ]` Implement the production lamination-parameter sensitivity path outside Abaqus native gradient support.
- `[ ]` Add restart/resume from checkpoint into the active driver flow.
- `[ ]` Add richer production logging/output formats.
- `[ ]` Add failure summaries and diagnostics tailored to real Abaqus runs.
- `[ ]` Revisit whether the newer dynamic solver stack should replace or wrap the legacy numerical core after the production path is stable.

## Current Snapshot

- `[x]` New pipeline/orchestration module exists.
- `[x]` Legacy numerical demos are now backed by regression-style tests.
- `[x]` Abaqus backend skeleton exists and is tested.
- `[x]` Full test suite currently passes.
- `[~]` The project now has the architectural scaffolding for the planned system, but the new pipeline still uses a generic gradient-penalty subproblem solver rather than the legacy min-max solver and laminate section model as its production backend.
- `[~]` The architecture now allows Abaqus-native gradients to be optional, but the production lamination-parameter sensitivity path still needs to be implemented explicitly on our side.
