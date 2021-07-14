/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.pistoncylinder;

import etomica.atom.IAtomKinetic;
import etomica.data.DataSourceScalar;
import etomica.integrator.IntegratorHardFasterer;
import etomica.space.Boundary;
import etomica.space.Vector;
import etomica.units.dimensions.Pressure;
import etomica.units.dimensions.Pressure2D;

/**
 * data source front for virial sum from P1HardMovingBoundary
 * returns pressure exerted on the wall by atoms
 */
public class DataSourceWallPressureFasterer extends DataSourceScalar implements IntegratorHardFasterer.CollisionListener {

    protected final IntegratorHardFasterer integrator;
    protected double collisionRadius;
    protected double virialSum;
    protected double lastTime;

    public DataSourceWallPressureFasterer(IntegratorHardFasterer integrator, double collisionRadius) {
        super("Wall Pressure", integrator.getBox().getSpace().D() == 2 ? Pressure2D.DIMENSION : Pressure.DIMENSION);
        this.integrator = integrator;
        integrator.addCollisionListener(this);
        virialSum = 0;
        lastTime = 0;
        this.collisionRadius = collisionRadius;
    }

    public void setCollisionRadius(double cr) {
        collisionRadius = cr;
    }

    public double getDataAsScalar() {
        double currentTime = integrator.getCurrentTime();
        double value = virialSum / (currentTime - lastTime);
        lastTime = currentTime;
        virialSum = 0;
        return value;
    }

    @Override
    public void pairCollision(IAtomKinetic atom1, IAtomKinetic atom2, Vector rij, Vector dv, double virial, double tCollision) {

    }

    @Override
    public void fieldCollision(IAtomKinetic atom, Vector r, Vector deltaP, double tCollision) {
        int dim = r.getD();
        double dp = deltaP.getX(1);
        Boundary b = integrator.getBox().getBoundary();
        double A = 1;
        for (int i = 0; i < dim; i++) {
            if (i == 1) continue;
            double L = b.getBoxSize().getX(i);
            A *= (L - 2 * collisionRadius);
        }
        if (dim == 2) {
            dp = -dp;
        }
        if (dp > 0) {
            virialSum += dp / A;
        }
    }
}
