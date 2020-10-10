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
import etomica.box.storage.Tokens;
import etomica.box.storage.VectorStorage;
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
    protected final Vector translationVector, oldPosition;
    protected final IRandom random;
    protected int iAtom;
    protected double uOld;
    protected double uNew = Double.NaN;
    protected AtomSource atomSource;
    protected boolean fixOverlap;
    protected Space space;
    protected final VectorStorage positions;

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
        oldPosition = space.makeVector();
        setStepSizeMax(stepSizeMax);
        setStepSizeMin(0.0);
        setStepSize(stepSize);
        perParticleFrequency = true;
        setBox(box);
        positions = box.getAtomStorage(Tokens.POSITION);
    }

    public boolean doTrial() {
        IAtom atom = atomSource.getAtom();
        if (atom == null) return false;
        iAtom = atom.getLeafIndex();
        uOld = potentialMasterFasterer.computeOneOld(iAtom);
        if (uOld > 1e8) {
            throw new RuntimeException("atom " + atom + " in box " + box + " has an overlap");
        }
        translationVector.setRandomCube(random);
        translationVector.TE(stepSize);
        Vector r = positions.get(iAtom);
        oldPosition.E(r);
        r.PE(translationVector);
        Vector shift = box.getBoundary().centralImage(r);
        r.PE(shift);
        potentialMasterFasterer.updateAtom(iAtom);
        return true;
    }//end of doTrial

    public double getChi(double temperature) {
        uNew = potentialMasterFasterer.computeOne(iAtom);
        return Math.exp(-(uNew - uOld) / temperature);
    }

    public double energyChange() {
        return uNew - uOld;
    }

    public void acceptNotify() {
//        System.out.println("accepted");
        potentialMasterFasterer.processAtomU(1);
        // put it back, then compute old contributions to energy
        Vector r = positions.get(iAtom);
        r.E(oldPosition);
        potentialMasterFasterer.updateAtom(iAtom);
        potentialMasterFasterer.computeOne(iAtom);
        r.PE(translationVector);
        Vector shift = box.getBoundary().centralImage(r);
        r.PE(shift);
        potentialMasterFasterer.processAtomU(-1);
        potentialMasterFasterer.updateAtom(iAtom);
    }

    public void rejectNotify() {
//        System.out.println("rejected");
        positions.get(iAtom).E(oldPosition);
        potentialMasterFasterer.updateAtom(iAtom);
    }

    public AtomIterator affectedAtoms() {
        affectedAtomIterator.setAtom(box.getLeafList().get(iAtom));
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
