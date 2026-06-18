/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.colloid;

import etomica.atom.AtomLeafAgentManager;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.potential.P2HardGeneric;
import etomica.species.ISpecies;

public class P2HardSphereMC extends P2HardGeneric {

    protected int chainLength;

    protected final AtomLeafAgentManager<? extends IAtomList> bondManager;
    protected double bondFac;
    protected final ISpecies speciesMonomer;

    public P2HardSphereMC(double sigma, double bondFac, int chainLength, AtomLeafAgentManager<? extends IAtomList> bondManager, ISpecies speciesMonomer) {
        super(new double[]{sigma}, new double[]{Double.POSITIVE_INFINITY}, true);
        this.bondManager = bondManager;
        this.speciesMonomer = speciesMonomer;
        setBondFac(bondFac);
        setChainLength(chainLength);
    }

    public void setChainLength(int newChainLength) {
        chainLength = newChainLength;
    }

    public double u(double r2) {
        return 0;
    }

    @Override
    protected double[] getEnergies(IAtom atom1, IAtom atom2) {
        int idxMonomer = atom1.getParentGroup().getType() == speciesMonomer ?
                atom1.getParentGroup().getIndex() :
                atom2.getParentGroup().getIndex();
        boolean grafted = idxMonomer % chainLength == 0;
        if (grafted) {
            return new double[]{Double.POSITIVE_INFINITY, 0, Double.POSITIVE_INFINITY};
        }
        return collisionDistances2;
    }

    @Override
    protected double[] getCollisionDistances2(IAtom atom1, IAtom atom2) {
        int idxMonomer = atom1.getParentGroup().getType() == speciesMonomer ?
                atom1.getParentGroup().getIndex() :
                atom2.getParentGroup().getIndex();
        boolean grafted = idxMonomer % chainLength == 0;
        if (grafted) {
            return new double[]{collisionDistances2[0] * bondFac * bondFac, collisionDistances2[0], collisionDistances2[0] + 1};
        }
        return collisionDistances2;
    }

    public void setBondFac(double newBondFac) {
        bondFac = newBondFac;
    }

    public double getBondFac() {
        return bondFac;
    }
}
