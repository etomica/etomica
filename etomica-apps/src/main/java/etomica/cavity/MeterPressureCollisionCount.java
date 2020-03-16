/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.cavity;

import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.integrator.IntegratorHard;
import etomica.units.dimensions.Pressure;

/**
 * Meter for the pressure by counting the number of collisions (assumes HS potential)
 */
public class MeterPressureCollisionCount extends DataSourceScalar implements IntegratorHard.CollisionListener {

    protected final int dim;
    protected long collisionCount = 0;
    protected final IntegratorHard integratorHard;
    protected double lastTime;

    public MeterPressureCollisionCount(IntegratorHard integrator) {
        super("Pressure(CC)", Pressure.dimension(integrator.getBox().getSpace().D()));
        integratorHard = integrator;
        integratorHard.addCollisionListener(this);
        lastTime = integratorHard.getCurrentTime();
        dim = integratorHard.getBox().getSpace().D();
    }

    public void reset() {
        collisionCount = 0;
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
        double value = (numAtomTemp + collisionCount * Math.sqrt(Math.PI * integratorHard.getTemperature()) / (dim * elapsedTime)) /
                box.getBoundary().volume();

        collisionCount = 0;
        lastTime = currentTime;
        return value;
    }

    /**
     * Implementation of CollisionListener interface
     * Adds collision virial (from potential) to accumulator
     */
    public void collisionAction(IntegratorHard.Agent agent) {
        collisionCount++;
    }

    public IntegratorHard getIntegrator() {
        return integratorHard;
    }
}
