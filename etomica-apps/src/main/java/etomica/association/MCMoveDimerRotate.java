/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.association;

import etomica.api.IAtom;
import etomica.api.IBox;
import etomica.api.IPotentialAtomic;
import etomica.api.IRandom;
import etomica.api.ISimulation;
import etomica.api.IVectorMutable;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomSource;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorArrayListSimple;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.mcmove.MCMoveBoxStep;
import etomica.nbr.cell.Api1ACell;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.space.ISpace;
import etomica.space.IVectorRandom;
import etomica.space.RotationTensor;

/**
 * Standard Monte Carlo atom-displacement trial move.
 *
 * @author David Kofke
 */
public class MCMoveDimerRotate extends MCMoveBoxStep {
    
    private static final long serialVersionUID = 2L;
    protected final AtomIteratorArrayListSimple affectedAtomIterator;
    protected final AtomArrayList affectedAtoms;
    protected final MeterPotentialEnergy energyMeter;
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
    protected transient RotationTensor rotationTensor;
    protected final IVectorMutable r0;
    protected AssociationManager associationManager;

    public MCMoveDimerRotate(ISimulation sim, PotentialMasterCell potentialMaster, ISpace _space, IPotentialAtomic dimerPotential) {
        this(potentialMaster, sim.getRandom(), _space, 1.0, 15.0, false, dimerPotential);
    }
    
    public MCMoveDimerRotate(PotentialMasterCell potentialMaster, IRandom random,
    		          ISpace _space, double stepSize, double stepSizeMax,
            boolean fixOverlap, IPotentialAtomic dimerPotential) {
        super(potentialMaster);
        this.affectedAtoms = new AtomArrayList(2);
        this.affectedAtomIterator = new AtomIteratorArrayListSimple(affectedAtoms);
        this.potentialMaster = potentialMaster;
        this.neighborIterator = new Api1ACell(3,1.0,potentialMaster.getCellAgentManager());
        this.dr = _space.makeVector();
        rotationTensor = _space.makeRotationTensor();
        this.r0= _space.makeVector();
        this.random = random;
        this.space = _space;
        this.dimerPotential = dimerPotential;
        atomSource = new AtomSourceRandomDimer();
        ((AtomSourceRandomDimer)atomSource).setRandomNumberGenerator(random);
        energyMeter = new MeterPotentialEnergy(potentialMaster);
        setStepSizeMax(stepSizeMax);
        setStepSizeMin(0.0);
        setStepSize(stepSize);
        perParticleFrequency = true;
        energyMeter.setIncludeLrc(false);
        this.fixOverlap = fixOverlap;
    }
    
    public void setAssociationManager(AssociationManager associationManager){
    	this.associationManager = associationManager;
    	((AtomSourceRandomDimer)atomSource).setAssociationManager(associationManager);
    }
    
    /**
     * Method to perform trial move.
     */
    public boolean doTrial() {
        atom = atomSource.getAtom();
        if (atom == null) return false;
        IAtom atom0 = atom;
        atom1 = associationManager.getAssociatedAtoms(atom).getAtom(0);
        energyMeter.setTarget(atom);
        uOld = energyMeter.getDataAsScalar();
        energyMeter.setTarget(atom1);
        uOld += energyMeter.getDataAsScalar();
        if(uOld > 1e8 && !fixOverlap) {
            throw new RuntimeException("atom "+atom+" in box "+box+" has an overlap");
        }
        double dTheta = (2*random.nextDouble() - 1.0)*stepSize;
        rotationTensor.setAxial(r0.getD() == 3 ? random.nextInt(3) : 2,dTheta);
        r0.E(atom1.getPosition());
        r0.PEa1Tv1(0.5,dr);
        doTransform(atom0);
        doTransform(atom1);
//        if (atom.getParentGroup().getIndex() == 371 || atom.getParentGroup().getIndex() == 224 ||
//        		((IAtom)atom1).getParentGroup().getIndex() == 371 || ((IAtom)atom1).getParentGroup().getIndex() == 224){
//        	System.out.println("MCMoveDimerRotate "+atom);
//        }
        
        return true;
    }//end of doTrial
    protected void doTransform(IAtom a) {
            IVectorMutable r = a.getPosition();
            r.ME(r0);
            box.getBoundary().nearestImage(r);
            rotationTensor.transform(r);
            r.PE(r0);
        }
    
    /**
     * Returns log of the ratio of the trial probabilities, ln(Tij/Tji) for the
     * states encountered before (i) and after (j) the most recent call to doTrial(). 
     * Tij is the probability that this move would generate state j from state i, and
     * Tji is the probability that a subsequent call to doTrial would return to state i
     * from state j.
     */
    public double getA() {
    	if (associationManager.getAssociatedAtoms(atom).getAtomCount() > 1) {
        	return 0;
        } 
        if (associationManager.getAssociatedAtoms(atom).getAtomCount() == 1){
        	IAtom atomj = associationManager.getAssociatedAtoms(atom).getAtom(0);
        	if(associationManager.getAssociatedAtoms(atomj).getAtomCount() > 1){
        		return 0;
        	} 
        }
    	return 1.0;
	}
    
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
        rotationTensor.invert();
        doTransform(atom);
        doTransform(atom1);
        
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
