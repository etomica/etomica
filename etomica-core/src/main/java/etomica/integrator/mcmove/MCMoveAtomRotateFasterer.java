/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.integrator.mcmove;

import etomica.atom.AtomSource;
import etomica.atom.AtomSourceRandomLeaf;
import etomica.atom.IAtomOriented;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorSinglet;
import etomica.box.Box;
import etomica.potential.compute.PotentialCompute;
import etomica.space.IOrientation;
import etomica.space.Space;
import etomica.space3d.Orientation3D;
import etomica.util.random.IRandom;

/**
 * Standard Monte Carlo atom-displacement trial move. Selects an atom at random, displaces it to a new position selected at random
 * on a cube centered on its current position. The set of atoms subject to the trial can be configured;
 * default is all leaf-atoms.
 *
 * @author David Kofke
 */
public class MCMoveAtomRotateFasterer extends MCMoveBoxStep {

    protected final PotentialCompute potentialCompute;
    protected final AtomIteratorSinglet affectedAtomIterator = new AtomIteratorSinglet();
    protected final IRandom random;
    protected IAtomOriented atom;
    protected double uOld;
    protected double uNew = Double.NaN;
    protected AtomSource atomSource;
    protected boolean fixOverlap;
    protected Space space;
    protected IOrientation oldOrientation, newOrientation;

    /**
     * Constructs the move with default stepSize = 1.0, stepSizeMax = 15.0, fixOverlap = false
     *
     * @param random           random number generator used to select the atom and its displacement
     * @param potentialCompute used to construct MeterPotentialEnergy required by full constructor
     */
    public MCMoveAtomRotateFasterer(IRandom random, PotentialCompute potentialCompute, Box box) {
        super();
        this.potentialCompute = potentialCompute;
        this.random = random;
        this.space = box.getSpace();
        atomSource = new AtomSourceRandomLeaf();
        ((AtomSourceRandomLeaf) atomSource).setRandomNumberGenerator(random);
        setStepSizeMax(Math.PI);
        setStepSizeMin(0.0);
        setStepSize(stepSize);
        perParticleFrequency = true;
        setBox(box);

        IAtomOriented atom0 = (IAtomOriented) box.getLeafList().get(0);

        if (atom0.getOrientation() instanceof Orientation3D) {
            oldOrientation = new Orientation3D(space);
            newOrientation = new Orientation3D(space);
        } else {
            oldOrientation = space.makeOrientation();
            newOrientation = space.makeOrientation();
        }

    }

    public boolean doTrial() {
        atom = (IAtomOriented) atomSource.getAtom();
        if (atom == null) return false;
        uOld = potentialCompute.computeOneOld(atom);
        if (uOld > 1e8) {
            throw new RuntimeException("atom " + atom + " in box " + box + " has an overlap");
        }

        IOrientation iOrientation = atom.getOrientation();
        oldOrientation.E(iOrientation);  //save old orientation
        iOrientation.randomRotation(random, stepSize);

        potentialCompute.updateAtom(atom);
        return true;
    }//end of doTrial

    public double getChi(double temperature) {
        uNew = potentialCompute.computeOne(atom);
        return Math.exp(-(uNew - uOld) / temperature);
    }

    public double energyChange() {
        return uNew - uOld;
    }

    public void acceptNotify() {
//        System.out.println("accepted");
        potentialCompute.processAtomU(1);

        // put it back, then compute old contributions to energy
        newOrientation.E(atom.getOrientation());
        atom.getOrientation().E(oldOrientation);

        potentialCompute.updateAtom(atom);
        potentialCompute.computeOne(atom);

        // put atom in new orientation again
        atom.getOrientation().E(newOrientation);

        potentialCompute.processAtomU(-1);
        potentialCompute.updateAtom(atom);
    }

    public void rejectNotify() {
//        System.out.println("rejected");
        atom.getOrientation().E(oldOrientation);
        potentialCompute.updateAtom(atom);
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
