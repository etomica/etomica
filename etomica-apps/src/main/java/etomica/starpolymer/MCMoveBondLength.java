/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.starpolymer;

import etomica.atom.AtomArrayList;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.integrator.mcmove.MCMoveMolecule;
import etomica.potential.compute.PotentialCompute;
import etomica.space.Vector;
import etomica.species.ISpecies;
import etomica.util.random.IRandom;

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
public class MCMoveBondLength extends MCMoveMolecule {

    protected Vector r0;
    protected double dr;
    protected int startAtom;
    protected ISpecies species;
    protected final int armLength;
    protected final Vector translateVector;
    protected double bl;
    protected boolean rejectme;
    protected final AtomArrayList movedAtoms;
    protected IAtom[] movedAtomArray;

    public MCMoveBondLength(PotentialCompute potentialMaster, IRandom random, Box box, int armLength) {
        this(potentialMaster, random, box, 1.0, armLength);
    }

    public MCMoveBondLength(PotentialCompute potentialMaster,
                            IRandom random, Box box, double stepSize, int armLength) {
        super(random, potentialMaster, box);
        this.armLength = armLength;
        translateVector = space.makeVector();
        movedAtoms = new AtomArrayList();
    }

    public void setSpecies(ISpecies newSpecies) {
        species = newSpecies;
    }

    //note that total energy is calculated
    public boolean doTrial() {
        molecule = moleculeSource.getMolecule();
        IAtomList childList = molecule.getChildList();
        int numChildren = childList.size();

        startAtom = random.nextInt(numChildren - 1) + 1;
        movedAtoms.clear();
        int arm = (startAtom-1) / armLength;
        for (int i=startAtom; i<1+(arm+1)*armLength; i++) {
            movedAtoms.add(childList.get(i));
        }

        IAtom fixedAtom = childList.get(startAtom - 1);
        if (startAtom % armLength == 1 || armLength == 1) {
//            System.out.println("star atom " + startAtom);
            fixedAtom = childList.get(0);
        }
        movedAtomArray = movedAtoms.toArray();
        uOld = potentialCompute.computeManyAtomsOld(movedAtomArray);
        if (uOld > 1e8) {
            throw new RuntimeException("molecule " + molecule + " in box " + box + " has an overlap");
        }

        r0 = fixedAtom.getPosition();
//            System.out.println(selectedAtoms[i]+" "+j+" before "+selectedAtoms[i].coord.position());
        dr = stepSize * 2 * (random.nextDouble() - 0.5);
        rejectme = false;
        doTransform();
        if (rejectme) return false;

        return true;
    }

    protected void doTransform() {
        IAtomList childList = molecule.getChildList();
        int stopAtom = (((startAtom - 1) / armLength) + 1) * armLength;
        Vector startPos = childList.get(startAtom).getPosition();
        translateVector.Ev1Mv2(startPos, r0);
        bl = Math.sqrt(translateVector.squared());
        if (dr < -bl) {
            rejectme = true;
            return;
        }
        translateVector.TE(dr / bl);
        for (int iChild = startAtom; iChild <= stopAtom; iChild++) {
            IAtom a = childList.get(iChild);
            Vector r = a.getPosition();
            r.PE(translateVector);
            r.PE(box.getBoundary().centralImage(r));
            potentialCompute.updateAtom(a);
        }
    }

    public void acceptNotify() {
        potentialCompute.processAtomU(1);
        // put it back, then compute old contributions to energy
        dr = -dr;
//        System.out.println("rejecting bond length move");
        doTransform();
        potentialCompute.computeManyAtoms(movedAtomArray);
        potentialCompute.processAtomU(-1);
        dr = -dr;
//        System.out.println("rejecting bond length move");
        doTransform();
    }

    public void rejectNotify() {
        dr = -dr;
//        System.out.println("rejecting bond length move");
        doTransform();
    }

    public double getChi(double temperature) {
        uNew = potentialCompute.computeManyAtoms(movedAtomArray);
        double blNew = bl + dr;
        double ratio = blNew / bl;
        return Math.exp(-(uNew - uOld) / temperature) * ratio * ratio;
    }
}
