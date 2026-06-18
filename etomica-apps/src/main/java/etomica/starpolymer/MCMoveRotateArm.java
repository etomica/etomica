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
import etomica.space3d.RotationTensor3D;
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
public class MCMoveRotateArm extends MCMoveMolecule {

    protected Vector r0;
    protected int startAtom;
    protected final Vector work1, work2, work3;
    protected ISpecies species;
    protected final int armLength;
    protected final RotationTensor3D rotationTensor;
    protected final Vector axis;
    protected final AtomArrayList movedAtoms;
    protected IAtom[] movedAtomArray;

    public MCMoveRotateArm(PotentialCompute potentialMaster, IRandom random, Box box, int armLength) {
        this(potentialMaster, random, box, 1.0, armLength);
    }

    public MCMoveRotateArm(PotentialCompute potentialMaster,
                           IRandom random, Box box, double stepSize, int armLength) {
        super(random, potentialMaster, box);
        this.armLength = armLength;
        setStepSizeMax(Math.PI);
        setStepSize(stepSize);
        work1 = space.makeVector();
        work2 = space.makeVector();
        work3 = space.makeVector();
        rotationTensor = new RotationTensor3D();
        axis = space.makeVector();
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
        movedAtomArray = movedAtoms.toArray();
        uOld = potentialCompute.computeManyAtomsOld(movedAtomArray);
        if (uOld > 1e8) {
            throw new RuntimeException("molecule " + molecule + " in box " + box + " has an overlap");
        }

        IAtom fixedAtom = childList.get(startAtom - 1);
        if (startAtom % armLength == 1 || armLength == 1) fixedAtom = childList.get(0);
        r0 = fixedAtom.getPosition();
//            System.out.println(selectedAtoms[i]+" "+j+" before "+selectedAtoms[i].coord.position());
        axis.setRandomSphere(random);
        double theta = random.nextDouble() * stepSize;
//        System.out.println("rotate from "+startAtom+" by "+theta+" about "+axis);
        rotationTensor.setRotationAxis(axis, theta);
        doTransform();

        return true;
    }


    protected void doTransform() {
        IAtomList childList = molecule.getChildList();
        int stopAtom = (((startAtom - 1) / armLength) + 1) * armLength;
        for (int iChild = startAtom; iChild <= stopAtom; iChild++) {
            IAtom a = childList.get(iChild);
            Vector r = a.getPosition();
            r.ME(r0);
            box.getBoundary().nearestImage(r);
            rotationTensor.transform(r);
            r.PE(r0);
            r.PE(box.getBoundary().centralImage(r));
            potentialCompute.updateAtom(a);
        }
    }

    public void acceptNotify() {
        potentialCompute.processAtomU(1);
        // put it back, then compute old contributions to energy
        rotationTensor.invert();
        doTransform();

        potentialCompute.computeManyAtoms(movedAtomArray);
        potentialCompute.processAtomU(-1);
        rotationTensor.invert();
        doTransform();
    }

    public void rejectNotify() {
//        System.out.println("rejecting rotate");
        rotationTensor.invert();
        doTransform();
    }

    public double getChi(double temperature) {
        uNew = potentialCompute.computeManyAtoms(movedAtomArray);

        return Math.exp(-(uNew - uOld) / temperature);
    }
}
