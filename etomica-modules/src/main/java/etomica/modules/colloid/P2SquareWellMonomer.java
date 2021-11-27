/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.colloid;

import etomica.atom.IAtom;
import etomica.atom.IAtomKinetic;
import etomica.atom.IAtomList;
import etomica.potential.IPotential2;
import etomica.potential.P2HardGeneric;
import etomica.space.Vector;

public class P2SquareWellMonomer extends P2HardGeneric implements IPotential2 {

    protected int chainLength;
    protected double rGraftMin2;

    public P2SquareWellMonomer(double sigma, double lambda, double epsMM, double bondFac, int chainLength) {
        super(new double[]{sigma, sigma * lambda}, new double[]{Double.POSITIVE_INFINITY, -epsMM}, true);
        setBondFac(bondFac);
        setChainLength(chainLength);
    }

    public void setRGraftMin(double rGraftMin) {
        rGraftMin2 = rGraftMin * rGraftMin;
    }

    public void setChainLength(int newChainLength) {
        chainLength = newChainLength;
    }

    public double u(Vector dr12, IAtom atom1, IAtom atom2) {
        int idx1 = atom1.getParentGroup().getIndex();
        int idx2 = atom2.getParentGroup().getIndex();
        int chainIdx1 = idx1 / chainLength;
        int chainIdx2 = idx2 / chainLength;
        int childIdx1 = idx1 % chainLength;
        int childIdx2 = idx2 % chainLength;
        boolean grafted = childIdx1 + childIdx2 == 0;
        boolean bonded = (chainIdx1 == chainIdx2) && Math.abs(childIdx1 - childIdx2) == 1;
        if (grafted || bonded) return 0;
        return u(dr12.squared());
    }

    public double energy(IAtomList pair) {
        throw new RuntimeException("nope");
    }

    @Override
    protected double[] getEnergies(IAtom atom1, IAtom atom2) {
        int idx1 = atom1.getParentGroup().getIndex();
        int idx2 = atom2.getParentGroup().getIndex();
        int childIdx1 = idx1 % chainLength;
        int childIdx2 = idx2 % chainLength;
        boolean grafted = childIdx1 + childIdx2 == 0;
        if (grafted) return new double[]{Double.POSITIVE_INFINITY};
        int chainIdx1 = idx1 / chainLength;
        int chainIdx2 = idx2 / chainLength;
        boolean bonded = (chainIdx1 == chainIdx2) && Math.abs(childIdx1 - childIdx2) == 1;
        if (bonded) return new double[]{Double.POSITIVE_INFINITY, 0, Double.POSITIVE_INFINITY};
        return energies;
    }

    @Override
    protected double[] getCollisionDistances2(IAtom atom1, IAtom atom2) {
        int idx1 = atom1.getParentGroup().getIndex();
        int idx2 = atom2.getParentGroup().getIndex();
        int childIdx1 = idx1 % chainLength;
        int childIdx2 = idx2 % chainLength;
        boolean grafted = childIdx1 + childIdx2 == 0;
        if (grafted) return new double[]{rGraftMin2};
        int chainIdx1 = idx1 / chainLength;
        int chainIdx2 = idx2 / chainLength;
        boolean bonded = (chainIdx1 == chainIdx2) && Math.abs(childIdx1 - childIdx2) == 1;
        if (bonded)
            return new double[]{bondFac * bondFac * collisionDistances2[0], collisionDistances2[0], collisionDistances2[0] + 1};
        return collisionDistances2;
    }

    protected int decideBump(IAtomKinetic atom1, IAtomKinetic atom2, int oldState, boolean core, double ke, double reducedMass, double bij, double r2, double[] du, double[] virial, double falseTime) {

        int idx1 = atom1.getParentGroup().getIndex();
        int idx2 = atom2.getParentGroup().getIndex();
        int chainIdx1 = idx1 / chainLength;
        int chainIdx2 = idx2 / chainLength;
        int childIdx1 = idx1 % chainLength;
        int childIdx2 = idx2 % chainLength;
        boolean grafted = childIdx1 + childIdx2 == 0;
        boolean bonded = (chainIdx1 == chainIdx2) && Math.abs(childIdx1 - childIdx2) == 1;
        if (grafted || bonded) {
            virial[0] = 2.0 * reducedMass * bij;
            return oldState;
        }
        return super.decideBump(atom1, atom2, oldState, core, ke, reducedMass, bij, r2, du, virial, falseTime);
    }

    public void setBondFac(double newBondFac) {
        bondFac = newBondFac;
    }

    public double getBondFac() {
        return bondFac;
    }

    protected double bondFac;
}
