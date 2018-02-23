/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data.meter;

import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.integrator.IntegratorHard;
import etomica.units.dimensions.Pressure;

/**
 * Meter for the pressure (given as the compressibility factor) of a hard potential.
 * Performs sum of collision virial over all collisions, and manipulates value
 * to obtain the compressibility factor, PV/NkT.
 *
 * @author David Kofke
 */
public class MeterPressureHard extends DataSourceScalar implements IntegratorHard.CollisionListener {

    protected final int dim;
    protected double virialSum = 0;
    protected final IntegratorHard integratorHard;
    protected double lastTime;

    public MeterPressureHard(IntegratorHard integrator) {
        super("Pressure", Pressure.dimension(integrator.getBox().getSpace().D()));
        integratorHard = integrator;
        integratorHard.addCollisionListener(this);
        lastTime = integratorHard.getCurrentTime();
        dim = integratorHard.getBox().getSpace().D();
    }

    public void reset() {
        virialSum = 0.0;
        lastTime = integratorHard.getCurrentTime();
    }

    /**
     * Returns P = (NT - (virial sum)/((elapsed time)*T*(space dimension)))/V
     * Virial sum and elapsed time apply to period since last call to this method.
     */
    public double getDataAsScalar() {
        Box box = integratorHard.getBox();
        double currentTime = integratorHard.getCurrentTime();
        double elapsedTime = currentTime - lastTime;
        if (elapsedTime == 0.0) return Double.NaN;
        if (elapsedTime < 0) throw new RuntimeException("you should have called reset");
        double numAtomTemp = integratorHard.getKineticEnergy() * 2 / dim;
        if (integratorHard.isIsothermal()) {
            numAtomTemp = integratorHard.getTemperature() * box.getLeafList().size();
        }
        double value = (numAtomTemp - virialSum / (dim * elapsedTime)) /
                box.getBoundary().volume();

        virialSum = 0.0;
        lastTime = currentTime;
        return value;
    }

    /**
     * Implementation of CollisionListener interface
     * Adds collision virial (from potential) to accumulator
     */
    public void collisionAction(IntegratorHard.Agent agent) {
        virialSum += agent.collisionPotential.lastCollisionVirial();
    }

    /**
     * Implementation of Meter.MeterCollisional interface.  Returns -(collision virial).
     * Suitable for tabulation of PV
     */
    public double collisionValue(IntegratorHard.Agent agent) {
        return -agent.collisionPotential.lastCollisionVirial();
    }

    public IntegratorHard getIntegrator() {
        return integratorHard;
    }
}
