# etomica-core quirks and suspected bugs

A running list of surprising behavior in etomica-core (and libraries it wraps), found while building
or changing code elsewhere. Entries are for David and co-developers to review and decide on. Nothing
here has been fixed unless the status says so. Maintained by Claude Code sessions (see CLAUDE.md).

Each entry gives where it is, what happens, how it was found, the workaround used (if any), a
possible fix, and a status: **open**, **decided: won't fix**, **fixed** (with commit), or
**unconfirmed** (seen once, not reproduced).

---

## 1. `DisplayPlotXChart` throws on paint when every series is empty

- **Where:** `etomica.graphics.DisplayPlotXChart` (`dataChanged`, around the `xMin > xMax` block
  that resets the axis limits to NaN).
- **What happens:** if all series have no points (for example a history cleared by Reinitialize
  while paused, or a plot shown before Start), the axis range is NaN and XChart throws
  `NumberFormatException` on every repaint. The plot stays blank and the console fills with traces.
- **Found:** Heat Capacity module, 2026-11.
- **Workaround:** the module seeds each history with one finite point after clearing it.
- **Possible fix:** when there are no finite points, keep the previous axis range (or a default such
  as 0–1) instead of NaN.
- **Status:** open. May affect other modules that clear histories while paused (not checked).

## 2. `DeviceDelaySlider` overrides any preset sleep period

- **Where:** `etomica.graphics.DeviceDelaySlider` constructor, which calls
  `delaySlider.setValue(0)`, and that sets the controller's sleep period to 0.
- **What happens:** a module that calls `controller.setSleepPeriod(n)` before building the GUI
  loses that setting once `SimulationGraphic` creates its delay slider. The slider also isn't
  exposed, so its value can't be set afterwards.
- **Found:** Heat Capacity module, 2026-11. The Thermal Conductivity module's
  `setSleepPeriod(1)` is probably also overridden (inferred from the code, not tested).
- **Workaround:** the Heat Capacity module makes its own delay slider.
- **Possible fix:** initialize the slider from `controller.getSleepPeriod()` (inverting the
  delay-exponent mapping) rather than to 0, and/or expose the slider.
- **Status:** open.

## 3. `DeviceSlider.setPrecision` can push a stale value through the modifier

- **Where:** `etomica.graphics.DeviceSlider.setPrecision`, which calls `doUpdate()` when a
  modifier is attached.
- **What happens:** with a modifier already attached (as `DeviceThermoSlider` has), changing the
  precision can apply a value that was rescaled for the old precision. In the Heat Capacity module
  this once set the integrator temperature to 0.01 at startup. It's timing-dependent, so it doesn't
  happen every run.
- **Found:** Heat Capacity module, 2026-11.
- **Workaround:** after `setPrecision`, restore T and attach a fresh
  `ModifierGeneral(integrator, "temperature")`.
- **Possible fix:** detach the modifier during `setPrecision` (as `setMinimum`/`setMaximum` already
  do), then resync the slider from the modifier's current value.
- **Status:** open.

## 4. `ANDERSEN_SINGLE` thermostat ignores `thermostatNoDrift`

- **Where:** `etomica.integrator.IntegratorMD.doThermostat`, `ANDERSEN_SINGLE` branch.
- **What happens:** the other thermostat branches call `shiftMomenta()` when `thermostatNoDrift`
  is set; `ANDERSEN_SINGLE` doesn't, so the net momentum random-walks. Removing that momentum
  afterwards (with `ActionZeroMomentum`) leaves the sampled temperature low by a fraction of about
  1/N_atoms: each randomization adds 3 dof of energy, and the zeroing then removes some of it.
- **Found:** Heat Capacity module, 2026-11 (measured average T about 1% low with 100 diatomics).
- **Workaround:** the module overrides `doThermostat` to redraw the relative velocity of a random
  atom pair, which conserves momentum exactly and samples the canonical distribution on the
  zero-momentum manifold (measured average T = 1.000).
- **Possible fix:** either document that `ANDERSEN_SINGLE` doesn't conserve momentum, or offer a
  momentum-conserving pair variant like the one above.
- **Status:** open.

## 5. XChart freezes the event thread when tick labels can't be made distinct

- **Where:** the XChart 3.6.6 axis tick calculator (library, not etomica), reached through
  `DisplayPlotXChart` when a decimal pattern is set on an axis.
- **What happens:** with an integer pattern such as `"0"` on an axis spanning less than about one
  unit, XChart keeps searching for distinct tick labels and the Swing event thread hangs.
- **Found:** Heat Capacity module, 2026-11 (time axis early in a run).
- **Workaround:** use a pattern with optional decimals, such as `"0.###"`.
- **Possible fix:** none in core is needed; could be noted in `DisplayPlotXChart` javadoc. Upgrading
  XChart is constrained (pinned to 3.6.6 for API reasons).
- **Status:** open (library behavior).

## 6. `DeviceSlider` tick labels have an uneven gap for some ranges

- **Where:** `etomica.graphics.DeviceSlider.setTicks`.
- **What happens:** for the Heat Capacity bond-stiffness slider (range 25–1000, 4 major ticks) the
  labels were unevenly spaced, leaving a visible gap. Cosmetic.
- **Found:** Heat Capacity module, 2026-11. Not investigated further.
- **Status:** open, cosmetic.

## 7. `BasicSliderUI` "highLabel is null" NullPointerException at startup

- **Where:** Swing `BasicSliderUI` painting a `DeviceSlider` with custom label tables; exact trigger
  unknown.
- **What happens:** one NPE trace at startup, after which the GUI worked normally.
- **Found:** Heat Capacity module, 2026-11. Seen once, not reproduced in later runs.
- **Status:** unconfirmed.

## 8. `NeighborIteratorCellFaster` applies periodic images to pairs involving a core atom

- **Where:** `etomica.nbr.cell.NeighborIteratorCellFaster` (added on the `Achiya` branch), the
  `nocelltype` (core) branches of `iterUpNeighbors`, `iterDownNeighbors` and `iterAllNeighbors`,
  which call `box.getBoundary().nearestImage(rij)`. Pairs between two non-core atoms use plain
  distances.
- **What happens:** `BoxCluster(cluster, space, L)` with L > 0 makes a periodic boundary. Mayer
  sampling needs isolated molecules, but every pair involving a core then uses its nearest-image
  distance once the two atoms are more than L/2 apart along any axis. The
  core potentials (`P2ParticleParticle`, `P2ParticleSegment`) have infinite range, so this always
  changes the energy, by up to about 1e-6 in u12 and 0.4 in the intramolecular energy in the
  checks below.
- **Related, outside core:** `MCMoveClusterRotateMoleculeMulti.doTransform` (etomica-apps) applies
  `nearestImage` to each atom's offset from the rotation centre and then `centralImage` to the new
  position, and `CenterOfMass.position` (core, correct for periodic systems) images too. In a
  periodic box this changes positions, not just energies: any atom more than L/2 from the centre,
  or a molecule straddling the box face, gets folded by a box length. In the periodic 60-monomer
  run, saved configurations had "bonds" up to 62 long, and the average intramolecular energy was
  2.3e8 against 179 with a non-periodic box. `MeterRadiusGyrationComponentsCore` (used by
  `VirialStarSingle`) also images its distances.
- **Found:** `VirialjackCellTarget` energy checks, 8 arms of 60 monomers, L = 50, rc = 2.5,
  2026-11-24, in a first version of that driver that made its box periodic.
- **Who is affected:** nothing found so far. `VirialjackCell` gets its boxes from
  `SimulationVirialOverlap2`, which makes them non-periodic (box lengths default to 0, and
  `VirialjackCell` only sets the size afterwards), so all the imaging calls are no-ops there. Any
  driver that does make a periodic box for these molecules would be affected.
- **Workaround:** `VirialjackCellTarget` uses a non-periodic boundary of the same size (`new
  BoxCluster(cluster, space)` then `setBoxSize`); nearest-image then does nothing. The branch's
  `NeighborCellManager` handles a non-periodic box (padding cells stay empty; atoms outside the box
  still go to the outer cells), and the energy checks then pass. This is the same set-up
  `VirialjackCell` already has.
- **Possible fix:** use a non-periodic boundary in the drivers, which fixes all of the above at
  once; or remove the imaging from each of these places.
- **Status:** open (latent; no affected driver known).

---

## Not quirks (for reference)

- **`NeighborListManager` "drifted into unsafe zone" warning:** this is correct behavior. In the
  Heat Capacity module it showed that a time step of 0.01 was too large for WCA collisions at
  T = 3, which the module fixed by lowering its time-step cap.
