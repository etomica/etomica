/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.pistoncylinder;

import etomica.data.meter.MeterPressureHard;
import etomica.integrator.IntegratorHard;
import etomica.potential.P1HardMovingBoundary;

/**
 * data source front for virial sum from P1HardMovingBoundary
 * returns pressure exerted on the wall by atoms
 */
public class DataSourceWallPressure extends MeterPressureHard {

    protected final P1HardMovingBoundary wallPotential;

    public DataSourceWallPressure(P1HardMovingBoundary pistonPotential, IntegratorHard integrator) {
        super(integrator);
        wallPotential = pistonPotential;
    }

    /**
     * Implementation of CollisionListener interface
     * Adds collision virial (from potential) to accumulator
     */
    public void collisionAction(IntegratorHard.Agent agent) {
        if (agent.collisionPotential == wallPotential) {
            if (dim == 2) {
                virialSum += wallPotential.lastWallVirial();
            } else {
                virialSum -= wallPotential.lastWallVirial();
            }
        }
    }

    public double getDataAsScalar() {
        double currentTime = integratorHard.getCurrentTime();
        double value = virialSum / (currentTime - lastTime);
        lastTime = currentTime;
        virialSum = 0;
        return value;
    }
}
