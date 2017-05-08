/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.freeenergy.npath;

import etomica.api.*;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomSource;
import etomica.atom.AtomSourceRandomLeaf;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorArrayListSimple;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.mcmove.MCMoveBoxStep;
import etomica.space.ISpace;
import etomica.space.IVectorRandom;

/**
 * Created by andrew on 4/11/17.
 */
public class MCMoveAtomCoupled extends MCMoveBoxStep {
    protected final AtomIteratorArrayListSimple affectedAtomIterator;
    protected final MeterPotentialEnergy energyMeter;
    protected final IVectorRandom translationVector;
    protected IAtom atom, atom2;
    protected double uOld;
    protected double uNew = Double.NaN;
    protected AtomSource atomSource;
    protected boolean fixOverlap;
    protected final IRandom random;
    protected ISpace space;
    protected IAtomList atoms;

    public MCMoveAtomCoupled(IRandom random, IPotentialMaster potentialMaster, ISpace _space) {
        this(potentialMaster, random, _space, 1.0, 15.0, false);
    }

    public MCMoveAtomCoupled(IPotentialMaster potentialMaster, IRandom random,
                      ISpace _space, double stepSize, double stepSizeMax,
                      boolean fixOverlap) {
        this(potentialMaster, new MeterPotentialEnergy(potentialMaster), random, _space,
                stepSize, stepSizeMax, fixOverlap);
    }

    public MCMoveAtomCoupled(IPotentialMaster potentialMaster, MeterPotentialEnergy meterPE, IRandom random,
                      ISpace _space, double stepSize, double stepSizeMax,
                      boolean fixOverlap) {
        super(potentialMaster);
        this.random = random;
        this.space = _space;
        atomSource = new AtomSourceRandomLeaf();
        ((AtomSourceRandomLeaf)atomSource).setRandomNumberGenerator(random);
        this.energyMeter = meterPE;
        translationVector = (IVectorRandom)space.makeVector();
        setStepSizeMax(stepSizeMax);
        setStepSizeMin(0.0);
        setStepSize(stepSize);
        perParticleFrequency = true;
        energyMeter.setIncludeLrc(false);
        this.fixOverlap = fixOverlap;
        affectedAtomIterator = new AtomIteratorArrayListSimple();
    }

    /**
     * Method to perform trial move.
     */
    public boolean doTrial() {
        atom = atomSource.getAtom();
        if (atom == null) return false;
        int n = atoms.getAtomCount();
        int idx2 = atom.getLeafIndex()+n/2;
        if (idx2>=n) idx2-=n;
        atom2 = atoms.getAtom(idx2);

        energyMeter.setTarget(atom);
        uOld = energyMeter.getDataAsScalar();
        energyMeter.setTarget(atom2);
        uOld += energyMeter.getDataAsScalar();
        if(uOld > 1e8 && !fixOverlap) {
            throw new RuntimeException("atom "+atom+" in box "+box+" has an overlap");
        }
        translationVector.setRandomCube(random);
        translationVector.TE(stepSize);
        atom.getPosition().PE(translationVector);
        atom2.getPosition().PE(translationVector);
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
        energyMeter.setTarget(atom);
        uNew = energyMeter.getDataAsScalar();
        energyMeter.setTarget(atom2);
        uNew += energyMeter.getDataAsScalar();
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
        atom2.getPosition().PE(translationVector);
    }

    public AtomIterator affectedAtoms() {
        AtomArrayList list = (AtomArrayList)affectedAtomIterator.getList();
        list.clear();
        list.add(atom);
        list.add(atom2);
        return affectedAtomIterator;
    }

    public void setBox(IBox p) {
        super.setBox(p);
        energyMeter.setBox(p);
        atomSource.setBox(p);
        atoms = p.getLeafList();
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
