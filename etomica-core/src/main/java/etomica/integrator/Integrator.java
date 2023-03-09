/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.integrator;

import etomica.space.Vector;
import etomica.util.Statefull;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.Writer;

/**
 * Integrator implements the algorithm used to move the atoms around and
 * generate new configurations in one or more boxs. All integrator techniques,
 * such as molecular dynamics or Monte Carlo, are implemented via subclasses of
 * this Integrator class. The Integrator's activities are managed via the
 * actions of the governing Controller.
 *
 * @author David Kofke and Andrew Schultz
 */
public abstract class Integrator implements Statefull {

    protected final IntegratorEventManager eventManager;
    protected boolean initialized = false;
    protected int interval;
    protected long stepCount;
    private int iieCount;

    public Integrator() {
        setEventInterval(1);
        stepCount = 0;
        eventManager = new IntegratorEventManager(this);
    }

    /**
     * @return value of interval (number of doStep calls) between
     * firing of interval events by integrator.
     */
    public int getEventInterval() {
        return interval;
    }

    /**
     * Sets value of interval between successive firing of integrator interval events.
     *
     * @param interval
     */
    public void setEventInterval(int interval) {
        this.interval = interval;
        if (iieCount > interval) {
            iieCount = interval;
        }
    }

    /**
     * Performs the elementary integration step, such as a molecular dynamics
     * time step, or a Monte Carlo trial.
     */
    public final void doStep() {
        stepCount++;
        --iieCount;
        if (iieCount == 0) {
            eventManager.stepStarted();
        }

        doStepInternal();

        if (iieCount == 0) {
            eventManager.stepFinished();
            iieCount = interval;
        }
    }

    /**
     * Performs the elementary integration step, such as a molecular dynamics
     * time step, or a Monte Carlo trial.
     */
    protected abstract void doStepInternal();

    /**
     * Returns the number of steps performed by the integrator since it was
     * initialized.
     */
    public long getStepCount() {
        return stepCount;
    }

    /**
     * Defines the actions taken by the integrator to reset itself, such as
     * required if a perturbation is applied to the simulated box (e.g.,
     * addition or deletion of a molecule). Also invoked when the
     * <p>
     * This should be called by subclasses before they have performed their own
     * reset integrator is started or initialized.
     */
    public void reset() {
        if (!initialized) {
            setup();
            initialized = true;
        }
        eventManager.initialized();
    }

    /**
     * @return true if reset() was ever called.
     */
    public boolean isInitialized() {
        return initialized;
    }

    /**
     * This method resets the step counter.
     */
    public void resetStepCount() {
        stepCount = 0;
    }


    public IntegratorEventManager getEventManager() {
        return eventManager;
    }

    /**
     * Perform initialization.  Subclasses can override this method to set up
     * before integration begins.
     */
    protected void setup() {
        resetStepCount();
        iieCount = interval;
    }

    /**
     * Integrator agent that holds a force vector. Used to indicate that an atom
     * could be under the influence of a force.
     */
    public interface Forcible {
        Vector force();
    }

    /**
     * Integrator agent that holds a torque vector.
     */
    public interface Torquable {
        Vector torque();
    }

    public void saveState(Writer fw) throws IOException {
        fw.write("" + stepCount + " " + iieCount + "\n");
    }

    public void restoreState(BufferedReader br) throws IOException {
        String[] bits = br.readLine().split(" ");
        stepCount = Long.parseLong(bits[0]);
        iieCount = Integer.parseInt(bits[1]);
    }

    public void postRestore() {
        eventManager.initialized();
    }
}
