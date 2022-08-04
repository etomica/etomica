/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.starpolymer;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.integrator.mcmove.MCMoveBoxStep;
import etomica.molecule.IMolecule;
import etomica.molecule.MoleculeSource;
import etomica.molecule.MoleculeSourceRandomMolecule;
import etomica.potential.compute.PotentialCompute;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.ISpecies;
import etomica.util.random.IRandom;
import etomica.virial.BoxCluster;

/**
 * An MC Move for cluster simulations that "wiggles" a chain molecule.  If the
 * first or last atom in the chain is chosen, it is moved to a new position
 * with the same bond length as before, but perturbed by some angle from its
 * original position.  If an Atom in the middle of the chain is chosen, a
 * crankshaft move is performed that maintains its distances with its
 * neighbors.  If a middle Atom has a bond angle too close to 180 degrees
 * (such that rotation does nothing) the Atom is not moved at all.
 * In each doTrial, wiggle moves are attempted on all molecules in the box.
 *
 * @author Andrew Schultz
 */
public class MCMoveClusterBondLength extends MCMoveBoxStep {

    protected IMolecule molecule;
    protected MoleculeSource moleculeSource;
    protected IRandom random;
    protected final PotentialCompute potentialMaster;
    protected Vector r0;
    protected double dr;
    protected int startAtom;
    protected final Space space;
    protected double wOld, wNew, uOld, uNew;
    protected ISpecies species;
    protected final int armLength;
    protected final Vector translateVector;
    protected double bl;

    public MCMoveClusterBondLength(IRandom random, PotentialCompute potentialMaster, int armLength, Space _space) {
        this(potentialMaster, random, 1.0, armLength, _space);
    }

    public MCMoveClusterBondLength(PotentialCompute potentialMaster,
                                   IRandom random, double stepSize, int armLength, Space _space) {
        super();
        moleculeSource = new MoleculeSourceRandomMolecule();
        ((MoleculeSourceRandomMolecule) moleculeSource).setRandomNumberGenerator(random);
        this.space = _space;
        this.armLength = armLength;
        translateVector = _space.makeVector();
        this.potentialMaster = potentialMaster;
    }

    public void setSpecies(ISpecies newSpecies) {
        species = newSpecies;
    }

    //note that total energy is calculated
    public boolean doTrial() {
        molecule = moleculeSource.getMolecule();
        uOld = potentialMaster.computeOneMolecule(molecule);
        wOld = ((BoxCluster) box).getSampleCluster().value((BoxCluster) box);
        IAtomList childList = molecule.getChildList();
        int numChildren = childList.size();

        startAtom = random.nextInt(numChildren - 1) + 1;

        IAtom fixedAtom = childList.get(startAtom - 1);
        if (startAtom % armLength == 1) {
//            System.out.println("star atom " + startAtom);
            fixedAtom = childList.get(0);
        }
        r0 = fixedAtom.getPosition();
//            System.out.println(selectedAtoms[i]+" "+j+" before "+selectedAtoms[i].coord.position());
        dr = stepSize * 2 * (random.nextDouble() - 0.5);
        doTransform();

        ((BoxCluster) box).trialNotify();
        return true;
    }

    protected void doTransform() {
        IAtomList childList = molecule.getChildList();
        int stopAtom = (((startAtom - 1) / armLength) + 1) * armLength;
        Vector startPos = childList.get(startAtom).getPosition();
        translateVector.Ev1Mv2(startPos, r0);
        bl = Math.sqrt(translateVector.squared());
        translateVector.TE(dr / bl);
        double fac = -((double) (stopAtom - startAtom)) / childList.size();
        for (int iChild = 0; iChild < childList.size(); iChild++) {
            IAtom a = childList.get(iChild);
            Vector r = a.getPosition();
            boolean foo = startAtom < iChild && iChild <= stopAtom;
            r.PEa1Tv1(foo ? (fac + 1) : fac, translateVector);
        }
    }

    public void rejectNotify() {
        dr = -dr;
//        System.out.println("rejecting bond length move");
        doTransform();
    }

    public double getChi(double temperature) {
        uNew = potentialMaster.computeOneMolecule(molecule);
        wNew = ((BoxCluster) box).getSampleCluster().value((BoxCluster) box);
        double blNew = bl + dr;
        double ratio = blNew / bl;
        return (wOld == 0 ? 1 : wNew / wOld) * Math.exp(-(uNew - uOld) / temperature) * ratio * ratio;
    }

    @Override
    public void acceptNotify() {

    }

    @Override
    public double energyChange() {
        return 0;
    }
}
