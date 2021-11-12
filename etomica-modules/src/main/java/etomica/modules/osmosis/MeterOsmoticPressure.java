/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.osmosis;

import etomica.atom.IAtomKinetic;
import etomica.data.DataSourceScalar;
import etomica.integrator.IntegratorHard;
import etomica.potential.P1HardBoundary;
import etomica.space.Vector;
import etomica.units.dimensions.Pressure;

/**
 * Osmotic pressure meter that calculates the difference in
 */
public class MeterOsmoticPressure extends DataSourceScalar implements IntegratorHard.CollisionListener {

    private double virialSum;
    private double collisionRadius;
    private final P1HardBoundary[] boundaryPotentials;
    private final IntegratorHard integrator;
    private double lastTime;

    public MeterOsmoticPressure(P1HardBoundary[] boundaryPotentials, IntegratorHard integrator) {
        super("Osmotic Pressure", Pressure.DIMENSION);
        this.boundaryPotentials = boundaryPotentials;
        this.integrator = integrator;
        integrator.addCollisionListener(this);
        lastTime = integrator.getCurrentTime();
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

    public double getDataAsScalar() {
        double currentTime = integrator.getCurrentTime();
        double value = virialSum / (currentTime - lastTime);
        lastTime = currentTime;
        virialSum = 0;

        // calculate accessible "area"
        Vector dimensions = integrator.getBox().getBoundary().getBoxSize();
        double area = 1;
        for (int i = 1; i < dimensions.getD(); i++) {
            area *= (dimensions.getX(i) - 2 * collisionRadius);
        }
        return value / area;
    }

    @Override
    public void pairCollision(IAtomKinetic atom1, IAtomKinetic atom2, Vector rij, Vector dv, double virial, double tCollision) {
    }

    @Override
    public void fieldCollision(IAtomKinetic atom, Vector r, Vector deltaP, double tCollision) {
        double x = r.getX(0);
        if (Math.abs(x) < 2 * collisionRadius) return;
        virialSum -= deltaP.getX(0);
    }
}
