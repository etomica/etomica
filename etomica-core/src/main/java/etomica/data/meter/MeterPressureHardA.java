/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.data.meter;

import etomica.atom.AtomType;
import etomica.atom.IAtomKinetic;
import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.integrator.IntegratorHard;
import etomica.space.Vector;
import etomica.units.dimensions.Pressure;

public class MeterPressureHardA extends DataSourceScalar implements IntegratorHard.CollisionListener {

    protected final AtomType typeA;
    protected final int dim;
    protected double virialSumBA = 0, virialSumAA = 0;
    protected final IntegratorHard integratorHard;
    protected double lastTime;

    public MeterPressureHardA(IntegratorHard integrator, AtomType typeA) {
        super("Pressure", Pressure.dimension(integrator.getBox().getSpace().D()));
        this.typeA = typeA;
        integratorHard = integrator;
        integratorHard.addCollisionListener(this);
        lastTime = integratorHard.getCurrentTime();
        dim = integratorHard.getBox().getSpace().D();
    }

    public void reset() {
        virialSumBA = virialSumAA = 0.0;
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
        lastTime = currentTime;

        double rv = (virialSumBA + 0.5*virialSumAA)/(dim * elapsedTime)/box.getBoundary().volume();
        virialSumBA = virialSumAA = 0.0;
        return rv;
    }

    @Override
    public void pairCollision(IAtomKinetic atom1, IAtomKinetic atom2, Vector r12, Vector dv, double virial, double tCollision) {
        if (atom1.getType() == typeA) {
            if (atom2.getType() == typeA) {
                virialSumAA += virial;
            }
            else {
                virialSumBA += virial;
            }
        }
        else if (atom2.getType() == typeA) {
            virialSumBA += virial;
        }
    }

    public IntegratorHard getIntegrator() {
        return integratorHard;
    }

}
