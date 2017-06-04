/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.multiharmonic;

import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.api.IRandom;
import etomica.space.Vector;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.integrator.mcmove.MCMoveBox;
import etomica.potential.P1Harmonic;

/**
 * MCMove which moves all atoms.  New coordinates are taken from a Gaussian
 * distribution with width and center taken to match the P1Harmonic that
 * governs the atom's potential.  All moves are accepted.
 * 
 * @author Andrew Schultz
 */
public class MCMoveMultiHarmonic extends MCMoveBox {

    public MCMoveMultiHarmonic(P1Harmonic p1, IRandom random) {
        super(null);
        this.p1 = p1;
        iterator = new AtomIteratorLeafAtoms();
        this.random = random;
    }
    
    public void setBox(Box newBox) {
        super.setBox(newBox);
        iterator.setBox(box);
    }

    public AtomIterator affectedAtoms() {
        return iterator;
    }

    public double energyChange() {
        return uNew-uOld;
    }

    public void acceptNotify() {
    }

    public boolean doTrial() {
        uOld = uNew;
        double s = p1.getSpringConstant();
        Vector x0 = p1.getX0();
        IAtomList atoms = box.getLeafList();
        uNew = 0;
        double sqrtS = Math.sqrt(s);
        for (int i=0; i<atoms.getAtomCount(); i++) {
            Vector p = atoms.getAtom(i).getPosition();
            for (int j=0; j<p.getD(); j++) {
                double r = random.nextGaussian()/sqrtS;
                uNew += 0.5*r*r*s;
                p.setX(j, r);
            }
            p.PE(x0);
        }
        return true;
    }

    public double getA() {
        return 1;
    }

    public double getB() {
        return 0;
    }

    public void rejectNotify() {
        throw new RuntimeException("oops");
    }

    private static final long serialVersionUID = 1L;
    protected final P1Harmonic p1;
    protected final AtomIteratorLeafAtoms iterator;
    protected final IRandom random;
    protected double uOld, uNew;
}
