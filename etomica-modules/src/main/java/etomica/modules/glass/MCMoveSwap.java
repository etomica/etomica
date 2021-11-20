/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.modules.glass;

import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.integrator.mcmove.MCMoveBox;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.potential.compute.PotentialCompute;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.ISpecies;
import etomica.util.random.IRandom;

/**
 * MC Move class that swaps the position of two atoms of different species.
 */
public class MCMoveSwap extends MCMoveBox {

    protected final PotentialCompute potentialCompute;
    protected IAtom atomA, atomB;
    protected double uOld, uNew;
    protected final ISpecies speciesA, speciesB;
    protected final IRandom random;
    protected final Vector tmp;

    public MCMoveSwap(Space space, IRandom random, PotentialCompute potentialMaster, ISpecies speciesA, ISpecies speciesB) {
        super();
        this.potentialCompute = potentialMaster;
        this.random = random;
        this.speciesA = speciesA;
        this.speciesB = speciesB;
        tmp = space.makeVector();
    }

    @Override
    public void setBox(Box box) {
        super.setBox(box);
    }

    @Override
    public double energyChange() {
        return uNew - uOld;
    }

    @Override
    public boolean doTrial() {
        IMoleculeList listA = box.getMoleculeList(speciesA);
        IMoleculeList listB = box.getMoleculeList(speciesB);
        if (listA.size() * listB.size() == 0) return false;
        IMolecule mol = listA.get(random.nextInt(listA.size()));
        atomA = mol.getChildList().get(0);
        mol = listB.get(random.nextInt(listB.size()));
        atomB = mol.getChildList().get(0);

        uOld = potentialCompute.computeManyAtomsOld(atomA, atomB);
        // since we are swapping, we don't have to worry about the direct
        // interaction being double-counted.  it's unchanged in the new
        // configuration and also double-counted there.
        if (uOld > 1e8) {
            throw new RuntimeException("atom " + atomA + " or " + atomB + " in box " + box + " has an overlap");
        }
        tmp.E(atomA.getPosition());
        atomA.getPosition().E(atomB.getPosition());
        atomB.getPosition().E(tmp);
        potentialCompute.updateAtom(atomA);
        potentialCompute.updateAtom(atomB);
        return true;
    }

    @Override
    public double getChi(double temperature) {
        uNew = potentialCompute.computeManyAtoms(atomA, atomB);
        return Math.exp(-(uNew - uOld) / temperature);
    }

    @Override
    public void acceptNotify() {
//        System.out.println("accepted");
        potentialCompute.processAtomU(1);
        // put it back, then compute old contributions to energy
        tmp.E(atomA.getPosition());
        atomA.getPosition().E(atomB.getPosition());
        atomB.getPosition().E(tmp);
        potentialCompute.updateAtom(atomA);
        potentialCompute.updateAtom(atomB);
        potentialCompute.computeManyAtoms(atomA, atomB);
        tmp.E(atomA.getPosition());
        atomA.getPosition().E(atomB.getPosition());
        atomB.getPosition().E(tmp);
        potentialCompute.processAtomU(-1);
        potentialCompute.updateAtom(atomA);
        potentialCompute.updateAtom(atomB);
    }

    @Override
    public void rejectNotify() {
        tmp.E(atomA.getPosition());
        atomA.getPosition().E(atomB.getPosition());
        atomB.getPosition().E(tmp);
        potentialCompute.updateAtom(atomA);
        potentialCompute.updateAtom(atomB);
    }
}
