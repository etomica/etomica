# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

Etomica is a molecular simulation framework in Java (Monte Carlo, molecular dynamics, virial coefficients), developed at the University at Buffalo. Javadocs: http://etomica.github.io/javadoc. The teaching modules in `etomica-modules` are described at http://modules.etomica.org (the site is sometimes down).

## Build & test

Gradle multi-project build (wrapper included), Java 17. All dependency and per-project config lives in the root `build.gradle`; there are no per-subproject build files except `benchmarks/build.gradle`.

Gradle 9 needs JDK 17+ to run, and `java_home -v 17` silently falls back to an older JDK if none is installed. On macOS, set `export JAVA_HOME=$(/usr/libexec/java_home -v 17)` before `./gradlew`. In IntelliJ, the project SDK, the compiler bytecode level and the Gradle JVM must all be 17, or builds fail with errors like "invalid target release".

```bash
./gradlew build                                   # compile + test everything (what CI runs, on Linux/macOS/Windows)
./gradlew :etomica-core:test                      # tests for one subproject
./gradlew :etomica-core:test --tests 'etomica.potential.EwaldSummationTest'          # one test class
./gradlew :etomica-core:test --tests 'etomica.nbr.*'                                 # a package
./gradlew :etomica-core:test --tests 'SomeTest.someMethod'                           # one method
```

Tests use JUnit 5 (Jupiter) and live in `etomica-core/src/test` and `etomica-apps/src/test`. Some tests touch Swing graphics; CI runs under xvfb.

Running a simulation class: any subproject gets the `application` plugin when `-Pmain` is given:

```bash
./gradlew :etomica-apps:run -Pmain=etomica.virial.simulations.VirialAlkane --args="-nPoints 3 -temperature 300"
./gradlew :etomica-core:run -Pmain=etomica.simulation.prototypes.LJMC3D
```

Other tasks:
- `./gradlew :etomica-modules:run` — launches the GUI teaching-module launcher (`etomica.modules.ModuleLauncher`); `:etomica-modules:archiveModulesApp` builds a jpackage'd distributable.
- `./benchmark.sh` — builds the JMH shadow jar (`:benchmarks:shadowJar`) and runs benchmarks; must be run from repo root.
- `./gradlew javadocAll` — aggregate javadoc into `docs/javadoc`.
- `-PbuildPython` — opt-in; creates a Python virtualenv in `etomica-core/src/main/resources/virtualenv` from `config/requirements.txt` (parmed_json) during `processResources`.

## Subprojects

- **etomica-core** — the framework: spaces/vectors, boxes, species, potentials, neighbor lists, integrators, MC moves, data pipeline, Swing graphics (`etomica.graphics`), and small prototype simulations (`etomica.simulation.prototypes`). Depends on etomica-graphics3D.
- **etomica-apps** — research simulations built on core, organized by project (`etomica.virial` is the largest — Mayer sampling/virial coefficients; also `normalmode`, `overlap`, `mappedvirial`, `models`, etc.). Many classes here are one-off research drivers.
- **etomica-modules** — interactive GUI teaching modules, registered with `ModuleLauncher`.
- **etomica-graph** — cluster-diagram (graph) algebra library used by the virial code.
- **etomica-graphics3D** — vendored 3D rendering (g3dsys + Jmol pieces).
- **benchmarks** — JMH benchmarks.
- `etomica-server`, `etomica_modules`, `stats_module*` directories are not part of the Gradle build.

## Core architecture

A simulation is a `Simulation` subclass that wires these pieces together in its constructor (see `etomica.simulation.prototypes.LJMC3D` as the canonical small example):

1. **Space** (`Space3D.getInstance()`, etc.) passed to `super(...)`; determines dimension and `Vector` implementation.
2. **Species** (`SpeciesGeneral`, usually built with `SpeciesBuilder` or `SpeciesGeneral.monatomic(...)`) added via `addSpecies` — must be added before boxes/potential masters that need the `SpeciesManager`. `AtomType`s belong to species.
3. **Box** (`addBox`), with a `Boundary`; populate with `box.setNMolecules`, size with `BoxInflate`, place with a `Configuration*` (e.g. `ConfigurationLattice`).
4. **Potentials / PotentialCompute** (`etomica.potential.compute.PotentialCompute`): the energy/force engine. Implementations include `PotentialMaster` (all pairs), `PotentialMasterCell` (cell lists), `PotentialMasterList` (neighbor lists, `etomica.nbr`), `PotentialMasterBonding` (intramolecular bonds via `BondingInfo`), `PotentialComputeField`, Ewald variants, and `PotentialComputeAggregate` to sum several. Pair potentials are registered per `AtomType` pair (`setPairPotential`); truncation is applied by wrapping (`P2SoftSphericalTruncated`, `TruncationFactory*`). Potentials report results through `PotentialCallback`s. Hard (discontinuous) potentials use the `*Hard` neighbor managers/iterators with `IntegratorHard`.
5. **Integrator** (`IntegratorMC` with `MCMove`s added to its `MCMoveManager`; or `IntegratorVelocityVerlet`, `IntegratorHard`, `IntegratorLangevin`, ...) constructed with the `PotentialCompute`.
6. **Controller**: `getController().addActivity(new ActivityIntegrate(integrator))`, or `runActivityBlocking(new ActivityIntegrate(integrator, steps))` for scripted runs (equilibrate, then reset and collect).

**Data pipeline** (`etomica.data`): `IDataSource` (meters in `etomica.data.meter`) → `DataPump`/`DataPumpListener` (added as an integrator listener with an interval) → `DataProcessor`s → `IDataSink`s (accumulators such as `AccumulatorAverageFixed`, `AccumulatorHistory`, or graphics displays). Data carries `IDataInfo` metadata; results from averaging accumulators are `DataGroup`s indexed by tags like `acc.AVERAGE.index`, `acc.ERROR.index`.

**Events**: integrators, boxes, and MC moves publish events through their own event managers (`IntegratorEventManager`, `BoxEventManager`, ...); listeners hook in there.

**Agent managers** (`AtomLeafAgentManager`, `BoxAgentManager`, `SpeciesAgentManager`) attach per-atom/box/species state and stay in sync as atoms are added/removed.

**Graphics**: `SimulationGraphic` wraps a simulation in a Swing UI; displays (`DisplayPlot`, `DisplayBox`, ...) act as data sinks. By convention, sims provide a nested `Graphic` class with its own `main` for the GUI version.

## Conventions

- Simulation parameters: a nested `public static class ... extends ParameterBase` with public fields; `main` calls `ParseArgs.doParseArgs(params, args)`, which maps `-fieldName value` (booleans: bare `-flag`; arrays: quoted `"0.3 0.5"`). When no args are given, `main` commonly sets parameters inline for interactive testing.
- Source files carry the MPL 2.0 header.
- `libs/ptolemy.jar` is a checked-in dependency (plotting).
- XChart is pinned to 3.6.6 deliberately (later versions change the API).

## Working with David

- I'm a co-lead developer and a domain expert (statistical mechanics, molecular simulation, Mayer sampling, mapped averaging).
- Active research areas: `etomica.virial` (Mayer-sampling MC for virial coefficients), mapped averaging (`etomica.mappedDensity`, `etomica.mappedvirial`, ...).
- Before any non-trivial change, propose a plan and wait for my approval.
- Prefer minimal, targeted edits. Don't refactor, rename, reformat, or "clean up" code outside the scope of the request.
- Match the style of the file being edited. Older research code in etomica-apps may use older patterns; don't modernize it unless asked.
- Prefer adding new classes alongside existing ones when changing behavior of a research driver.
- Don't change shared core classes (etomica-core) without flagging it first; many apps depend on them.
- Verification: after a change, compile the affected subproject.
- Never commit or push; I review and commit myself.

## Core quirks list

`CORE_QUIRKS.md` (repo root) records surprising behavior and suspected bugs in etomica-core (and
libraries it wraps) for me and my co-developer to review. Keep it current:

- When you hit such behavior during any task, especially when you work around it outside core, add an
  entry following the file's format: where, what happens, how it was found (module/task and date),
  the workaround, a possible fix, and a status.
- Before adding, check for an existing entry and update it instead (for example, more places it shows
  up, or a confirmed reproduction).
- Don't fix listed items in core unless I ask; the list is for deciding. When an item is fixed or I
  decide against fixing it, update its status rather than deleting it.
- Mention in your final report when you've added or changed entries.