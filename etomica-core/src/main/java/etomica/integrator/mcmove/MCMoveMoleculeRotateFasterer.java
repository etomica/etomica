/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.integrator.mcmove;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorSinglet;
import etomica.box.Box;
import etomica.molecule.CenterOfMass;
import etomica.molecule.IMolecule;
import etomica.molecule.MoleculeSource;
import etomica.molecule.MoleculeSourceRandomMolecule;
import etomica.potential.compute.PotentialCompute;
import etomica.space.RotationTensor;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.util.random.IRandom;

import java.util.ArrayList;
import java.util.List;

/**
 * Move that rotates a molecule about its center of mass.
 */
public class MCMoveMoleculeRotateFasterer extends MCMoveBoxStep {

    protected final PotentialCompute potentialCompute;
    protected final AtomIteratorSinglet affectedAtomIterator = new AtomIteratorSinglet();
    private final List<Vector> oldPositions = new ArrayList<>();
    protected final IRandom random;
    protected IMolecule molecule;
    protected double uOld;
    protected double uNew = Double.NaN;
    protected MoleculeSource moleculeSource;
    protected boolean fixOverlap;
    protected Space space;
    protected final Vector r0;
    protected final RotationTensor rotationTensor;

    /**
     * Constructs the move with default stepSize = 1.0, stepSizeMax = 15.0, fixOverlap = false
     *
     * @param random           random number generator used to select the atom and its displacement
     * @param potentialCompute used to construct MeterPotentialEnergy required by full constructor
     */
    public MCMoveMoleculeRotateFasterer(IRandom random, PotentialCompute potentialCompute, Box box) {
        super(null);
        this.potentialCompute = potentialCompute;
        this.random = random;
        this.space = box.getSpace();
        MoleculeSourceRandomMolecule source = new MoleculeSourceRandomMolecule();
        source.setRandomNumberGenerator(random);
        moleculeSource = source;
        setStepSizeMax(Math.PI);
        setStepSizeMin(0.0);
        setStepSize(Math.PI / 4);
        perParticleFrequency = true;
        setBox(box);
        r0 = space.makeVector();
        rotationTensor = space.makeRotationTensor();
    }

    public boolean doTrial() {
        molecule = moleculeSource.getMolecule();
        if (molecule == null) return false;
        uOld = potentialCompute.computeOneOldMolecule(molecule);
        if (uOld > 1e8) {
            throw new RuntimeException("molecule " + molecule + " in box " + box + " has an overlap");
        }
        while (oldPositions.size() < molecule.getChildList().size()) {
            oldPositions.add(space.makeVector());
        }
        for (IAtom a : molecule.getChildList()) {
            oldPositions.get(a.getIndex()).E(a.getPosition());
        }
        double dTheta = (2 * random.nextDouble() - 1.0) * stepSize;
        rotationTensor.setAxial(r0.getD() == 3 ? random.nextInt(3) : 2, dTheta);

        r0.E(CenterOfMass.position(box, molecule));
        doTransform();

        return true;
    }//end of doTrial

    protected void doTransform() {
        IAtomList childList = molecule.getChildList();
        for (IAtom a : childList) {
            Vector r = a.getPosition();
            r.ME(r0);
            box.getBoundary().nearestImage(r);
            rotationTensor.transform(r);
            r.PE(r0);
            r.PE(box.getBoundary().centralImage(r));
            potentialCompute.updateAtom(a);
        }
    }

    public double getChi(double temperature) {
        uNew = potentialCompute.computeOneMolecule(molecule);
//        System.out.println("rotate "+molecule.getIndex()+" "+uOld+" => "+uNew);
        return Math.exp(-(uNew - uOld) / temperature);
    }

    public double energyChange() {
        return uNew - uOld;
    }

    public void acceptNotify() {
//        System.out.println("accepted");
        potentialCompute.processAtomU(1);
        // put it back, then compute old contributions to energy
        molecule.getChildList().forEach(atom -> {
            atom.getPosition().E(oldPositions.get(atom.getIndex()));
            potentialCompute.updateAtom(atom);
        });
        potentialCompute.computeOneMolecule(molecule);
        potentialCompute.processAtomU(-1);
        doTransform();
    }

    public void rejectNotify() {
//        System.out.println("rejected");
        molecule.getChildList().forEach(atom -> {
            atom.getPosition().E(oldPositions.get(atom.getIndex()));
            potentialCompute.updateAtom(atom);
        });
    }

    public AtomIterator affectedAtoms() {
        throw new UnsupportedOperationException();
    }

    public void setBox(Box p) {
        super.setBox(p);
        moleculeSource.setBox(p);
    }

    /**
     * The MoleculeSource is used to select the molecule at the beginning of the trial
     *
     * @return the molecule source.
     */
    public MoleculeSource getMoleculeSource() {
        return moleculeSource;
    }

    /**
     * @param source The molecule source to set.
     */
    public void setMoleculeSource(MoleculeSource source) {
        moleculeSource = source;
    }
}
