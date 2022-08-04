/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.multiharmonic;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveBox;
import etomica.potential.P1Harmonic;
import etomica.space.Vector;
import etomica.util.random.IRandom;

/**
 * MCMove which moves all atoms.  New coordinates are taken from a Gaussian
 * distribution with width and center taken to match the P1Harmonic that
 * governs the atom's potential.  All moves are accepted.
 *
 * @author Andrew Schultz
 */
public class MCMoveMultiHarmonic extends MCMoveBox {

    public MCMoveMultiHarmonic(IntegratorMC integratorMC, P1Harmonic p1, IRandom random) {
        super();
        this.integratorMC = integratorMC;
        this.p1 = p1;
        this.random = random;
        uNew = Double.NaN;
    }

    public double energyChange() {
        return uNew-uOld;
    }

    public void acceptNotify() {
    }

    public boolean doTrial() {
        uOld = integratorMC.getPotentialEnergy();
        double s = p1.getSpringConstant();
        double sqrtS = Math.sqrt(s);
        Vector x0 = p1.getX0();
        IAtomList atoms = box.getLeafList();
        uNew = 0;
        for (IAtom atom : atoms) {
            Vector p = atom.getPosition();
            for (int j = 0; j < p.getD(); j++) {
                double r = random.nextGaussian() / sqrtS;
                uNew += 0.5 * r * r * s;
                p.setX(j, r);
            }
            p.PE(x0);
        }
        return true;
    }

    public double getChi(double temperature) {
        return 1;
    }

    public void rejectNotify() {
        throw new RuntimeException("oops");
    }

    protected final IntegratorMC integratorMC;
    protected final P1Harmonic p1;
    protected final IRandom random;
    protected double uOld, uNew;
}
