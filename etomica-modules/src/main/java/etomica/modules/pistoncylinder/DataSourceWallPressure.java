/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.pistoncylinder;

import etomica.data.meter.MeterPressureHard;
import etomica.integrator.IntegratorHard;
import etomica.potential.P1HardMovingBoundary;
import etomica.space.Space;

/**
 * data source front for virial sum from P1HardMovingBoundary
 * returns pressure exerted on the wall by atoms
 */
public class DataSourceWallPressure extends MeterPressureHard {

    public DataSourceWallPressure(Space _space, P1HardMovingBoundary pistonPotential) {
        super(_space);
        wallPotential = pistonPotential;
        space = _space;
    }
    
    /**
     * Implementation of CollisionListener interface
     * Adds collision virial (from potential) to accumulator
     */
    public void collisionAction(IntegratorHard.Agent agent) {
        if (agent.collisionPotential == wallPotential) {
            if (space.D() == 2) {
                virialSum += wallPotential.lastWallVirial();
            }
            else {
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
    
    private static final long serialVersionUID = 1L;
    protected final P1HardMovingBoundary wallPotential;
    private Space space;
}
