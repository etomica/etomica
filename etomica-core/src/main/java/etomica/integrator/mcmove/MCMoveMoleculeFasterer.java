/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.integrator.mcmove;

import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorSinglet;
import etomica.box.Box;
import etomica.molecule.IMolecule;
import etomica.molecule.MoleculeSource;
import etomica.molecule.MoleculeSourceRandomMolecule;
import etomica.potential.compute.PotentialCompute;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.util.random.IRandom;

import java.util.ArrayList;
import java.util.List;

/**
 * Standard Monte Carlo atom-displacement trial move. Selects an atom at random, displaces it to a new position selected at random
 * on a cube centered on its current position. The set of atoms subject to the trial can be configured;
 * default is all leaf-atoms.
 *
 * @author David Kofke
 */
public class MCMoveMoleculeFasterer extends MCMoveBoxStep {

    protected final PotentialCompute potentialCompute;
    protected final AtomIteratorSinglet affectedAtomIterator = new AtomIteratorSinglet();
    protected final Vector translationVector;
    private final List<Vector> oldPositions = new ArrayList<>();
    protected final IRandom random;
    protected IMolecule molecule;
    protected double uOld;
    protected double uNew = Double.NaN;
    protected MoleculeSource moleculeSource;
    protected boolean fixOverlap;
    protected Space space;

    /**
     * Constructs the move with default stepSize = 1.0, stepSizeMax = 15.0, fixOverlap = false
     *  @param random          random number generator used to select the atom and its displacement
     * @param potentialCompute used to construct MeterPotentialEnergy required by full constructor
     */
    public MCMoveMoleculeFasterer(IRandom random, PotentialCompute potentialCompute, Box box) {
        super(null);
        this.potentialCompute = potentialCompute;
        this.random = random;
        this.space = box.getSpace();
        MoleculeSourceRandomMolecule source = new MoleculeSourceRandomMolecule();
        source.setRandomNumberGenerator(random);
        moleculeSource = source;
        translationVector = space.makeVector();
        setStepSizeMax(box.getBoundary().getBoxSize().getX(0) / 2);
        setStepSizeMin(0.0);
        perParticleFrequency = true;
        setBox(box);
    }

    public boolean doTrial() {
        molecule = moleculeSource.getMolecule();
        if (molecule == null) return false;
        uOld = potentialCompute.computeOneOldMolecule(molecule);
        if (uOld > 1e8) {
            throw new RuntimeException("molecule " + molecule + " in box " + box + " has an overlap");
        }
        translationVector.setRandomCube(random);
        translationVector.TE(stepSize);
        while (oldPositions.size() < molecule.getChildList().size()) {
            oldPositions.add(space.makeVector());
        }
        molecule.getChildList().forEach(atom -> {
            oldPositions.get(atom.getIndex()).E(atom.getPosition());
            atom.getPosition().PE(translationVector);
            Vector shift = box.getBoundary().centralImage(atom.getPosition());
            atom.getPosition().PE(shift);
            potentialCompute.updateAtom(atom);
        });
        return true;
    }//end of doTrial

    public double getChi(double temperature) {
        uNew = potentialCompute.computeOneMolecule(molecule);
//        System.out.println("translate "+molecule.getIndex()+" "+uOld+" => "+uNew);
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
        molecule.getChildList().forEach(atom -> {
            atom.getPosition().PE(translationVector);
            Vector shift = box.getBoundary().centralImage(atom.getPosition());
            atom.getPosition().PE(shift);
            potentialCompute.updateAtom(atom);
        });
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
