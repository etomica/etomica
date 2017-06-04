/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.integrator.mcmove;

import etomica.atom.IAtom;
import etomica.potential.PotentialMaster;
import etomica.space.Vector;
import etomica.box.Box;
import etomica.api.IRandom;
import etomica.atom.AtomSource;
import etomica.atom.AtomSourceRandomLeaf;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorSinglet;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.space.Space;

/**
 * Standard Monte Carlo atom-displacement trial move.
 *
 * @author David Kofke
 */
public class MCMoveAtom extends MCMoveBoxStep {

    protected final AtomIteratorSinglet affectedAtomIterator = new AtomIteratorSinglet();
    protected final MeterPotentialEnergy energyMeter;
    protected final Vector translationVector;
    protected IAtom atom;
    protected double uOld;
    protected double uNew = Double.NaN;
    protected AtomSource atomSource;
    protected boolean fixOverlap;
    protected final IRandom random;
    protected Space space;

    public MCMoveAtom(IRandom random, PotentialMaster potentialMaster, Space _space) {
        this(potentialMaster, random, _space, 1.0, 15.0, false);
    }

    public MCMoveAtom(PotentialMaster potentialMaster, IRandom random,
                      Space _space, double stepSize, double stepSizeMax,
                      boolean fixOverlap) {
        this(potentialMaster, new MeterPotentialEnergy(potentialMaster), random, _space,
                stepSize, stepSizeMax, fixOverlap);
    }

    public MCMoveAtom(PotentialMaster potentialMaster, MeterPotentialEnergy meterPE, IRandom random,
                      Space _space, double stepSize, double stepSizeMax,
                      boolean fixOverlap) {
        super(potentialMaster);
        this.random = random;
        this.space = _space;
        atomSource = new AtomSourceRandomLeaf();
        ((AtomSourceRandomLeaf)atomSource).setRandomNumberGenerator(random);
        this.energyMeter = meterPE;
        translationVector = space.makeVector();
        setStepSizeMax(stepSizeMax);
        setStepSizeMin(0.0);
        setStepSize(stepSize);
        perParticleFrequency = true;
        energyMeter.setIncludeLrc(false);
        this.fixOverlap = fixOverlap;
    }

    /**
     * Method to perform trial move.
     */
    public boolean doTrial() {
        atom = atomSource.getAtom();
        if (atom == null) return false;
        energyMeter.setTarget(atom);
        uOld = energyMeter.getDataAsScalar();
        if(uOld > 1e8 && !fixOverlap) {
            throw new RuntimeException("atom "+atom+" in box "+box+" has an overlap");
        }
        translationVector.setRandomCube(random);
        translationVector.TE(stepSize);
        atom.getPosition().PE(translationVector);
        return true;
    }//end of doTrial

    /**
     * Returns log of the ratio of the trial probabilities, ln(Tij/Tji) for the
     * states encountered before (i) and after (j) the most recent call to doTrial(). 
     * Tij is the probability that this move would generate state j from state i, and
     * Tji is the probability that a subsequent call to doTrial would return to state i
     * from state j.
     */
    public double getA() {return 1.0;}

    /**
     * Returns the log of the limiting-distribution probabilities of states, ln(Pj/Pi), 
     * for the states encountered before (i) and after (j) the most recent call to 
     * doTrial.
     */
    public double getB() {
        uNew = energyMeter.getDataAsScalar();
        return -(uNew - uOld);
    }

    public double energyChange() {return uNew - uOld;}

    /**
     * Method called by IntegratorMC in the event that the most recent trial is accepted.
     */
    public void acceptNotify() {  /* do nothing */
    }

    /**
     * Method called by IntegratorMC in the event that the most recent trial move is
     * rejected.  This method should cause the system to be restored to the condition
     * before the most recent call to doTrial.
     */
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
     * @return Returns the atomSource.
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
