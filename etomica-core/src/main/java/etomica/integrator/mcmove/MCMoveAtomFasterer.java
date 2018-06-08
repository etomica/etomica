/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.integrator.mcmove;

import etomica.atom.AtomSource;
import etomica.atom.AtomSourceRandomLeaf;
import etomica.atom.IAtom;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorSinglet;
import etomica.box.Box;
import etomica.potential.PotentialMasterFasterer;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.util.random.IRandom;

/**
 * Standard Monte Carlo atom-displacement trial move. Selects an atom at random, displaces it to a new position selected at random
 * on a cube centered on its current position. The set of atoms subject to the trial can be configured;
 * default is all leaf-atoms.
 *
 * @author David Kofke
 */
public class MCMoveAtomFasterer extends MCMoveBoxStep {

    protected final PotentialMasterFasterer potentialMasterFasterer;
    protected final AtomIteratorSinglet affectedAtomIterator = new AtomIteratorSinglet();
    protected final Vector translationVector;
    protected final IRandom random;
    protected IAtom atom;
    protected double uOld;
    protected double uNew = Double.NaN;
    protected AtomSource atomSource;
    protected boolean fixOverlap;
    protected Space space;

    /**
     * Constructs the move with default stepSize = 1.0, stepSizeMax = 15.0, fixOverlap = false
     *
     * @param random          random number generator used to select the atom and its displacement
     * @param potentialMaster used to construct MeterPotentialEnergy required by full constructor
     */
    public MCMoveAtomFasterer(IRandom random, PotentialMasterFasterer potentialMaster, Box box) {
        super(null);
        this.potentialMasterFasterer = potentialMaster;
        this.random = random;
        this.space = box.getSpace();
        atomSource = new AtomSourceRandomLeaf();
        ((AtomSourceRandomLeaf) atomSource).setRandomNumberGenerator(random);
        translationVector = space.makeVector();
        setStepSizeMax(stepSizeMax);
        setStepSizeMin(0.0);
        setStepSize(stepSize);
        perParticleFrequency = true;
        setBox(box);
    }

    public boolean doTrial() {
        atom = atomSource.getAtom();
        if (atom == null) return false;
        uOld = potentialMasterFasterer.computeOneOld(atom);
        if (uOld > 1e8) {
            throw new RuntimeException("atom " + atom + " in box " + box + " has an overlap");
        }
        translationVector.setRandomCube(random);
        translationVector.TE(stepSize);
        atom.getPosition().PE(translationVector);
        return true;
    }//end of doTrial

    public double getChi(double temperature) {
        uNew = potentialMasterFasterer.computeOne(atom);
        return Math.exp(-(uNew - uOld) / temperature);
    }

    public double energyChange() {
        return uNew - uOld;
    }

    public void acceptNotify() {  /* do nothing */
//        System.out.println("accepted");
        potentialMasterFasterer.processAtomU(1);
        // put it back, then compute old contributions to energy
        translationVector.TE(-1);
        atom.getPosition().PE(translationVector);
        potentialMasterFasterer.computeOne(atom);
        translationVector.TE(-1);
        atom.getPosition().PE(translationVector);
        potentialMasterFasterer.processAtomU(-1);
    }

    public void rejectNotify() {
//        System.out.println("rejected");
        translationVector.TE(-1);
        atom.getPosition().PE(translationVector);
        potentialMasterFasterer.resetAtomDU();
    }

    public AtomIterator affectedAtoms() {
        affectedAtomIterator.setAtom(atom);
        return affectedAtomIterator;
    }

    public void setBox(Box p) {
        super.setBox(p);
        atomSource.setBox(p);
    }

    /**
     * The AtomSource is used to select the atom at the beginning of the trial
     *
     * @return the atomSource.
     */
    public AtomSource getAtomSource() {
        return atomSource;
    }

    /**
     * @param source The atomSource to set.
     */
    public void setAtomSource(AtomSource source) {
        atomSource = source;
    }
}
