/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.integrator.mcmove;

import etomica.box.Box;

/**
 * MCMove parent class for moves that operate on a single Box.
 *
 * @author Andrew Schultz
 */
public abstract class MCMoveBox extends MCMove {

    public MCMoveBox() {
        this(new MCMoveTracker());
    }

    public MCMoveBox(MCMoveTracker acceptanceTracker) {
        super(acceptanceTracker);
        perParticleFrequency = false;
    }

    /**
     * Sets the box on which this move acts.  The box itself can be changed
     * via this method, if desired.
     */
    public void setBox(Box p) {
        box = p;
    }

    /**
     * @return the box(s) on which this move acts.
     */
    public Box getBox() {
        return box;
    }

    public double energyChange(Box aBox) {
        if (box == aBox) {
            return energyChange();
        }
        return 0;
    }

    public abstract double energyChange();

    /**
     * Indicates whether this move should nominally be performed at a frequency
     * proportional to the number of molecules in the box.
     *
     * @see #nominalFrequency
     */
    public final boolean isNominallyPerParticleFrequency() {
        return perParticleFrequency;
    }

    public String toString() {
        return super.toString() +" "+box;
    }

    /**
     * Flag indicating whether nominal frequency is interpreted as a
     * perParticleFrequency, or as a full frequency. Default is false, but may
     * be given a different value by subclasses.
     */
    protected boolean perParticleFrequency;

    protected Box box;
}
