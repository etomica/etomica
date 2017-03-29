/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.integrator.mcmove;

import etomica.api.IAtom;
import etomica.api.IAtomList;
import etomica.api.IBox;
import etomica.api.IPotentialAtomic;
import etomica.api.IRandom;
import etomica.api.ISimulation;
import etomica.api.IVectorMutable;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomSource;
import etomica.atom.AtomSourceRandomLeaf;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorArrayListSimple;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.nbr.cell.Api1ACell;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.space.ISpace;
import etomica.space.IVectorRandom;

/**
 * Standard Monte Carlo atom-displacement trial move.
 *
 * @author David Kofke
 */
public class MCMoveDimer extends MCMoveBoxStep {
    
    private static final long serialVersionUID = 2L;
    protected final AtomIteratorArrayListSimple affectedAtomIterator;
    protected final AtomArrayList affectedAtoms;
    protected final MeterPotentialEnergy energyMeter;
    protected final IVectorRandom translationVector;
    protected IAtom atom;
    protected double uOld;
    protected double uNew = Double.NaN;
    protected AtomSource atomSource;
    protected boolean fixOverlap;
    protected final IRandom random;
    protected ISpace space;
    protected final PotentialMasterCell potentialMaster;
    protected final Api1ACell neighborIterator;
    protected final IVectorMutable dr;
    protected final IPotentialAtomic dimerPotential;
    protected IAtom atom1;

    public MCMoveDimer(ISimulation sim, PotentialMasterCell potentialMaster, ISpace _space, IPotentialAtomic dimerPotential) {
        this(potentialMaster, sim.getRandom(), _space, 1.0, 15.0, false, dimerPotential);
    }
    
    public MCMoveDimer(PotentialMasterCell potentialMaster, IRandom random,
    		          ISpace _space, double stepSize, double stepSizeMax,
            boolean fixOverlap, IPotentialAtomic dimerPotential) {
        super(potentialMaster);
        this.affectedAtoms = new AtomArrayList(2);
        this.affectedAtomIterator = new AtomIteratorArrayListSimple(affectedAtoms);
        this.potentialMaster = potentialMaster;
        this.neighborIterator = new Api1ACell(3,1.0,potentialMaster.getCellAgentManager());
        this.dr = _space.makeVector();
        this.random = random;
        this.space = _space;
        this.dimerPotential = dimerPotential;
        atomSource = new AtomSourceRandomLeaf();
        ((AtomSourceRandomLeaf)atomSource).setRandomNumberGenerator(random);
        energyMeter = new MeterPotentialEnergy(potentialMaster);
        translationVector = (IVectorRandom)space.makeVector();
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
        neighborIterator.setBox(box);
        neighborIterator.setTarget(atom);
        neighborIterator.reset();
        IAtom atom0 = atom;
        atom1 = null;
        for (IAtomList pair = neighborIterator.next(); pair != null;
             pair = neighborIterator.next()) {
        	 atom1 = pair.getAtom(1);
        	dr.Ev1Mv2(atom0.getPosition(),atom1.getPosition());
        	box.getBoundary().nearestImage(dr);
        	if (dr.squared() < 1) {
        		 double u = dimerPotential.energy(pair);
        		 if (u < 0) {
        			 break;
        		 }
        	}
        	atom1 = null;
        	        }
        if (atom1 == null){
        	return false;
        }
        energyMeter.setTarget(atom);
        uOld = energyMeter.getDataAsScalar();
        energyMeter.setTarget(atom1);
        uOld += energyMeter.getDataAsScalar();
        if(uOld > 1e8 && !fixOverlap) {
            throw new RuntimeException("atom "+atom+" in box "+box+" has an overlap");
        }
        translationVector.setRandomCube(random);
        translationVector.TE(stepSize);
        atom.getPosition().PE(translationVector);
        atom1.getPosition().PE(translationVector);
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
        energyMeter.setTarget(atom);
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
        atom1.getPosition().PE(translationVector);
    }
        
    
    public AtomIterator affectedAtoms() {
        affectedAtoms.clear();
        affectedAtoms.add(atom);
        affectedAtoms.add(atom1);
        return affectedAtomIterator;
    }
    
    public void setBox(IBox p) {
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
