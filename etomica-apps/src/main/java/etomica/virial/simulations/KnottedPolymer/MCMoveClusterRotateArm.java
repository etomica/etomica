/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations.KnottedPolymer;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.molecule.MoleculeSourceRandomMolecule;
import etomica.potential.PotentialCalculationEnergySum;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.RotationTensor3D;
import etomica.space3d.Vector3D;
import etomica.species.ISpecies;
import etomica.util.random.IRandom;
import etomica.virial.BoxCluster;
import etomica.virial.MCMoveClusterMolecule;

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
public class MCMoveClusterRotateArm extends MCMoveClusterMolecule {

    public PotentialMaster potentialMaster;

    public MCMoveClusterRotateArm(Simulation sim, PotentialMaster potentialMaster, int armLength, Space _space) {
        this(potentialMaster, sim.getRandom(), 1.0, armLength, _space);
    }

    public MCMoveClusterRotateArm(PotentialMaster potentialMaster,
                                  IRandom random, double stepSize, int armLength, Space _space) {
        super(random, _space, stepSize);
        this.potentialMaster = potentialMaster;
        moleculeSource = new MoleculeSourceRandomMolecule();
        ((MoleculeSourceRandomMolecule) moleculeSource).setRandomNumberGenerator(random);
        this.space = _space;
        this.armLength = armLength;
        setStepSizeMax(Math.PI);
        energyMeter = new MeterPotentialEnergy(potentialMaster);
        work1 = _space.makeVector();
        work2 = _space.makeVector();
        work3 = _space.makeVector();
        rotationTensor = new RotationTensor3D();
        axis = space.makeVector();
    }

    public void setBox(Box p) {
        super.setBox(p);
        energyMeter.setBox(p);
    }

    public void setSpecies(ISpecies newSpecies) {
        species = newSpecies;
    }

    //note that total energy is calculated
    public boolean doTrial() {
        molecule = moleculeSource.getMolecule();
        energyMeter.setTarget(molecule);
        uOld = energyMeter.getDataAsScalar();
        wOld = ((BoxCluster) box).getSampleCluster().value((BoxCluster) box);
        if (uOld > 1e8) {
            PotentialCalculationEnergySum.debug = true;
            energyMeter.getDataAsScalar();
            throw new RuntimeException("molecule " + molecule + " in box " + box + " has an overlap");
        }
        IAtomList childList = molecule.getChildList();
        int numChildren = childList.size();

        startAtom = random.nextInt(numChildren - 1) + 1;

        IAtom fixedAtom = childList.get(startAtom - 1);
        if (startAtom % armLength == 1) fixedAtom = childList.get(0);
        r0 = fixedAtom.getPosition();
//            System.out.println(selectedAtoms[i]+" "+j+" before "+selectedAtoms[i].coord.position());
        axis.setRandomSphere(random);
        double theta = random.nextDouble() * stepSize;
//        System.out.println("rotate from "+startAtom+" by "+theta+" about "+axis);
        rotationTensor.setRotationAxis(axis, theta);
        doTransform();
        ((BoxCluster) box).trialNotify();
        return true;
    }


    protected void doTransform() {
        IAtomList childList = molecule.getChildList();
        int stopAtom = (((startAtom - 1) / armLength) + 1) * armLength;
        Vector sum = new Vector3D();
        for (int iChild = startAtom; iChild <= stopAtom; iChild++) {
            IAtom a = childList.get(iChild);
            Vector r = a.getPosition();
            sum.ME(r);
            r.ME(r0);
            box.getBoundary().nearestImage(r);
            rotationTensor.transform(r);
            r.PE(r0);
            sum.PE(r);
        }
        sum.TE(1.0 / childList.size());
        for (int i = 0; i < childList.size(); i++) {
            IAtom a = childList.get(i);
            Vector r = a.getPosition();
            r.ME(sum);
        }
    }

    public void acceptNotify() {
//        System.out.println("accepting rotate");
        super.acceptNotify();
    }

    public void rejectNotify() {
//        System.out.println("rejecting rotate");
        rotationTensor.invert();
        doTransform();
        super.rejectNotify();
    }

    public double getChi(double temperature) {
        uNew = energyMeter.getDataAsScalar();
        wNew = ((BoxCluster) box).getSampleCluster().value((BoxCluster) box);
        return (wOld == 0 ? 1 : wNew / wOld) * Math.exp(-(uNew - uOld) / temperature);
    }

    protected final MeterPotentialEnergy energyMeter;
    protected Vector r0;
    protected int startAtom;
    protected final Vector work1, work2, work3;
    protected final Space space;
    protected ISpecies species;
    protected double wOld, wNew;
    protected final int armLength;
    protected final RotationTensor3D rotationTensor;
    protected final Vector axis;
}
