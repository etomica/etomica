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

---

## Not quirks (for reference)

- **`NeighborListManager` "drifted into unsafe zone" warning:** this is correct behavior. In the
  Heat Capacity module it showed that a time step of 0.01 was too large for WCA collisions at
  T = 3, which the module fixed by lowering its time-step cap.
