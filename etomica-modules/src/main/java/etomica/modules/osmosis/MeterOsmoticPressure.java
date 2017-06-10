/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.osmosis;
import etomica.space.Vector;
import etomica.data.meter.MeterPressureHard;
import etomica.integrator.IntegratorHard;
import etomica.potential.P1HardBoundary;
import etomica.space.Space;

/**
 * Osmotic pressure meter that calculates the difference in 
 */
public class MeterOsmoticPressure extends MeterPressureHard {
    
    private static final long serialVersionUID = 1L;
    private double virialSum;
    private double collisionRadius;
    private final P1HardBoundary[] boundaryPotentials;
    
    public MeterOsmoticPressure(Space space, P1HardBoundary[] boundaryPotentials) {
        super(space);
        this.boundaryPotentials = boundaryPotentials;
    }

    /**
     * Sets the collision radius used to calculate accessible "area"
     * assuming that the relevant hard "boundaries" form right-angles with
     * other hard boundaries.  The given collision radius should be the 
     * collsion radius of the molecules with the adjacent boundaries.  If
     * the adjacent boundaries are periodic, this should be 0.
     */
    public void setCollisionRadius(double newRadius) {
        collisionRadius = newRadius;
    }
    
    /**
     * Returns the collision radius used to calculate the accessible "area".
     */
    public double getCollisionRadius() {
        return collisionRadius;
    }
    
    public void collisionAction(IntegratorHard.Agent agent) {
        for (int i=0; i<boundaryPotentials.length; i++) {
            if (agent.collisionPotential == boundaryPotentials[i]) {
                virialSum += boundaryPotentials[i].lastCollisionVirialTensor().component(0, 0);
            }
        }
    }

    public double getDataAsScalar() {
        double currentTime = integratorHard.getCurrentTime();
        double value = virialSum / (currentTime - lastTime);
        lastTime = currentTime;
        virialSum = 0;

        // calculate accessible "area"
        Vector dimensions = integratorHard.getBox().getBoundary().getBoxSize();
        double area = 1;
        for (int i=1; i<dimensions.getD(); i++) {
            area *= (dimensions.getX(i)-2*collisionRadius);
        }
        return value / area;
    }
}
