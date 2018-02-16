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
import etomica.data.meter.MeterPotentialEnergy;
import etomica.potential.PotentialCalculationEnergySum;
import etomica.potential.PotentialMaster;
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
public class MCMoveAtom extends MCMoveBoxStep {

    protected final AtomIteratorSinglet affectedAtomIterator = new AtomIteratorSinglet();
    protected final MeterPotentialEnergy energyMeter;
    protected final Vector translationVector;
    protected final IRandom random;
    protected final boolean fixOverlap;
    protected final Space space;
    protected IAtom atom;
    protected double uOld;
    protected double uNew = Double.NaN;
    protected AtomSource atomSource;

    /**
     * Constructs the move with default stepSize = 1.0, stepSizeMax = 15.0, fixOverlap = false
     *
     * @param random          random number generator used to select the atom and its displacement
     * @param potentialMaster used to construct MeterPotentialEnergy required by full constructor
     * @param space           space of the simulation
     */
    public MCMoveAtom(IRandom random, PotentialMaster potentialMaster, Space space) {
        this(random, potentialMaster, space, 1.0, 15.0, false);
    }

    /**
     * @param random          random number generator used to select the atom and its displacement
     * @param potentialMaster used to construct MeterPotentialEnergy required by full constructor
     * @param space           space of the simulation
     * @param stepSize        starting step size for the trial
     * @param stepSizeMax     maximum allowable value of stepSize
     * @param fixOverlap      flag specifying whether trial throws an exception if the configuration
     *                        at the start of the trial has an overlap
     */
    public MCMoveAtom(IRandom random, PotentialMaster potentialMaster,
                      Space space, double stepSize, double stepSizeMax,
                      boolean fixOverlap) {
        this(random, potentialMaster, space, stepSize, stepSizeMax, fixOverlap,
                new MeterPotentialEnergy(potentialMaster));
    }

    /**
     * Constructs with specification of special-purpose energy meter.
     * @param random          random number generator used to select the atom and its displacement
     * @param potentialMaster not directly used
     * @param space           space of the simulation
     * @param stepSize        starting step size for the trial
     * @param stepSizeMax     maximum allowable value of stepSize
     * @param fixOverlap      if false, throws an exception if the configuration
     *                        at the start of the trial has an overlap
     * @param meterPE         used to compute energies before and after displacement
     */
    public MCMoveAtom(IRandom random, PotentialMaster potentialMaster,
                      Space space, double stepSize, double stepSizeMax,
                      boolean fixOverlap, MeterPotentialEnergy meterPE) {
        super(potentialMaster);
        this.random = random;
        this.space = space;
        atomSource = new AtomSourceRandomLeaf();
        ((AtomSourceRandomLeaf) atomSource).setRandomNumberGenerator(random);
        this.energyMeter = meterPE;
        translationVector = space.makeVector();
        setStepSizeMax(stepSizeMax);
        setStepSizeMin(0.0);
        setStepSize(stepSize);
        perParticleFrequency = true;
        energyMeter.setIncludeLrc(false);
        this.fixOverlap = fixOverlap;
    }

    public boolean doTrial() {
         atom = atomSource.getAtom();
        if (atom == null) return false;
        if(potential.isPotentialHard() && !fixOverlap) {
            uOld = 0.0;
        } else {
            energyMeter.setTarget(atom);
            uOld = energyMeter.getDataAsScalar();
            if (uOld > 1e8 && !fixOverlap) {
                PotentialCalculationEnergySum.debug = true;
                uOld = energyMeter.getDataAsScalar();
                throw new RuntimeException("atom " + atom + " in box " + box + " has an overlap");
            }
        }

        translationVector.setRandomCube(random);
        translationVector.TE(stepSize);
        atom.getPosition().PE(translationVector);
        return true;
    }//end of doTrial

    public double getChi(double temperature) {
        energyMeter.setTarget(atom);
        uNew = energyMeter.getDataAsScalar();
        return Math.exp(-(uNew - uOld) / temperature);
    }

    public double energyChange() {
        return uNew - uOld;
    }

    public void acceptNotify() {  /* do nothing */
    }

    public void rejectNotify() {
        translationVector.TE(-1);
        atom.getPosition().PE(translationVector);
    }

    public AtomIterator affectedAtoms() {
        affectedAtomIterator.setAtom(atom);
        return affectedAtomIterator;
    }

    public void setBox(Box p) {
        super.setBox(p);
        energyMeter.setBox(p);
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
