/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.association;

import etomica.atom.AtomArrayList;
import etomica.atom.AtomSource;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorArrayListSimple;
import etomica.box.Box;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.mcmove.MCMoveBoxStep;
import etomica.nbr.cell.Api1ACell;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.potential.IPotentialAtomic;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.util.random.IRandom;

/**
 * Standard Monte Carlo atom-displacement trial move of Trimer
 *
 * @author Hye Min Kim
 */
public class MCMoveSmer extends MCMoveBoxStep {
    
    private static final long serialVersionUID = 2L;
    protected final AtomIteratorArrayListSimple affectedAtomIterator;
    protected final MeterPotentialEnergy energyMeter;
    protected final Vector translationVector;
    protected IAtom atom;
    protected double uOld;
    protected double uNew = Double.NaN;
    protected AtomSource atomSource;
    protected boolean fixOverlap;
    protected final IRandom random;
    protected Space space;
    protected final PotentialMasterCell potentialMaster;
    protected final Api1ACell neighborIterator;
    protected final Vector dr;
    protected final IPotentialAtomic trimerPotential;
    protected final AtomArrayList smerList;
    protected final AtomArrayList newSmerList;
    protected AssociationManager associationManager;

    public MCMoveSmer(Simulation sim, PotentialMasterCell potentialMaster, Space _space, IPotentialAtomic dimerPotential) {
        this(potentialMaster, sim.getRandom(), _space, 1.0, 15.0, false, dimerPotential);
    }
    
    public MCMoveSmer(PotentialMasterCell potentialMaster, IRandom random,
                      Space _space, double stepSize, double stepSizeMax,
                      boolean fixOverlap, IPotentialAtomic TrimerPotential) {
        super(potentialMaster);
        this.smerList = new AtomArrayList();
        this.newSmerList = new AtomArrayList();
        this.affectedAtomIterator = new AtomIteratorArrayListSimple(smerList);
        this.potentialMaster = potentialMaster;
        this.neighborIterator = new Api1ACell(3,1.0,potentialMaster.getCellAgentManager());
        this.dr = _space.makeVector();
        this.random = random;
        this.space = _space;
        this.trimerPotential = TrimerPotential;
        atomSource = new AtomSourceRandomDimer();
        ((AtomSourceRandomDimer)atomSource).setRandomNumberGenerator(random);
        energyMeter = new MeterPotentialEnergy(potentialMaster);
        translationVector = space.makeVector();
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
        populateList(smerList);
//        if (smerList.indexOf(box.getLeafList().getAtom(388))>-1 ||smerList.indexOf(box.getLeafList().getAtom(115))>-1 ){
//        	System.out.println("moving smerList "+smerList);
//        	System.out.println("position 388 "+((IAtomPositioned)box.getLeafList().getAtom(388)).getPosition()+"position 115 "+((IAtomPositioned)box.getLeafList().getAtom(115)).getPosition());
//        }
        uOld = 0.0;
        for (int i=0; i<smerList.getAtomCount(); i+=1){
        	energyMeter.setTarget(smerList.getAtom(i));
            uOld += energyMeter.getDataAsScalar();
            if(uOld > 1e8 && !fixOverlap) {
            	//PotentialCalculationEnergySum.dodebug = true;
            	energyMeter.getDataAsScalar();
                throw new RuntimeException("smerList.getAtom(i) "+smerList.getAtom(i)+" in box "+box+" has an overlap");
            }
        }     
        if(uOld > 1e8 && !fixOverlap) {
            throw new RuntimeException("atom "+atom+" in box "+box+" has an overlap");
        }
        translationVector.setRandomCube(random);
        translationVector.TE(stepSize);
        for (int i=0; i<smerList.getAtomCount(); i+=1){
        	(smerList.getAtom(i)).getPosition().PE(translationVector);
        }
        
//        if (atom.getParentGroup().getIndex() == 371 || atom.getParentGroup().getIndex() == 224 ||
//        		((IAtom)atom1).getParentGroup().getIndex() == 371 || ((IAtom)atom1).getParentGroup().getIndex() == 224){
//        	System.out.println("MCMoveDimer "+atom);
//        }
        return true;
    }//end of doTrial
    
    protected int populateList(AtomArrayList mySmerList){
    	mySmerList.clear();
    	mySmerList.add(atom);
    	IAtomList bondList = associationManager.getAssociatedAtoms(atom);
    	if (bondList.getAtomCount() > 2){
    		return 0;
    	}
    	if (bondList.getAtomCount() == 2){
    		IAtom atom0 = bondList.getAtom(0);
    		IAtom atom1 = bondList.getAtom(1);
    		dr.Ev1Mv2((atom0).getPosition(), (atom1).getPosition());//dr = distance from the atom0 to atom1
        	box.getBoundary().nearestImage(dr);
        	double innerRadius = 0.8;
        	double minDistance = 2*(innerRadius*innerRadius)*(1+Math.cos(etomica.units.Degree.UNIT.toSim(27.0)));
        	if (dr.squared() < minDistance){
        		return 0;
        	}
    	}
    	if (bondList.getAtomCount() == 0){
    		return 1;
    	}
    	IAtom thisAtom = bondList.getAtom(0);
    	mySmerList.add(thisAtom);
    	IAtomList bondList1 = associationManager.getAssociatedAtoms(thisAtom);
    	if (bondList1.getAtomCount() > 2){
    		return 0;
    	}
    	if (bondList1.getAtomCount() == 2){
    		IAtom atom0 = bondList1.getAtom(0);
    		IAtom atom1 = bondList1.getAtom(1);
    		dr.Ev1Mv2((atom0).getPosition(), (atom1).getPosition());//dr = distance from the atom0 to atom1
        	box.getBoundary().nearestImage(dr);
        	double innerRadius = 0.8;
        	double minDistance = 2*(innerRadius*innerRadius)*(1+Math.cos(etomica.units.Degree.UNIT.toSim(27.0)));
        	if (dr.squared() < minDistance){
        		return 0;
        	}
    	}
    	IAtom previousAtom = atom;
    	while (bondList1.getAtomCount() > 1){
    		IAtom nextAtom = bondList1.getAtom(0);
    		if (nextAtom == previousAtom){
    			nextAtom = bondList1.getAtom(1);
    		} 
    		if (nextAtom == atom){
    			return 1;
    		}
    		mySmerList.add(nextAtom);
    		bondList1 = associationManager.getAssociatedAtoms(nextAtom);
    		if (bondList1.getAtomCount() > 2){
    			return 0;
        	}
    		if (bondList1.getAtomCount() == 2){
        		IAtom atom0 = bondList1.getAtom(0);
        		IAtom atom1 = bondList1.getAtom(1);
        		dr.Ev1Mv2((atom0).getPosition(), (atom1).getPosition());//dr = distance from the atom0 to atom1
            	box.getBoundary().nearestImage(dr);
            	double innerRadius = 0.8;
            	double minDistance = 2*(innerRadius*innerRadius)*(1+Math.cos(etomica.units.Degree.UNIT.toSim(27.0)));
            	if (dr.squared() < minDistance){
            		return 0;
            	}
        	}
    		previousAtom = thisAtom;
    		thisAtom = nextAtom;
    	}
    	if (bondList.getAtomCount()>1){
    		thisAtom = bondList.getAtom(1);
        	mySmerList.add(thisAtom);
        	bondList1 = associationManager.getAssociatedAtoms(thisAtom);
        	if (bondList1.getAtomCount() > 2){
    			return 0;
        	}
        	if (bondList1.getAtomCount() == 2){
        		IAtom atom0 = bondList1.getAtom(0);
        		IAtom atom1 = bondList1.getAtom(1);
        		dr.Ev1Mv2((atom0).getPosition(), (atom1).getPosition());//dr = distance from the atom0 to atom1
            	box.getBoundary().nearestImage(dr);
            	double innerRadius = 0.8;
            	double minDistance = 2*(innerRadius*innerRadius)*(1+Math.cos(etomica.units.Degree.UNIT.toSim(27.0)));
            	if (dr.squared() < minDistance){
            		return 0;
            	}
        	}
        	previousAtom = atom;
        	while (bondList1.getAtomCount() > 1){
        		IAtom nextAtom = bondList1.getAtom(0);
        		if (nextAtom == previousAtom){
        			nextAtom = bondList1.getAtom(1);
        		} 
        		mySmerList.add(nextAtom);
        		bondList1 = associationManager.getAssociatedAtoms(nextAtom);
        		if (bondList1.getAtomCount() > 2){
        			return 0;
            	}
        		if (bondList1.getAtomCount() == 2){
            		IAtom atom0 = bondList1.getAtom(0);
            		IAtom atom1 = bondList1.getAtom(1);
            		dr.Ev1Mv2((atom0).getPosition(), (atom1).getPosition());//dr = distance from the atom0 to atom1
                	box.getBoundary().nearestImage(dr);
                	double innerRadius = 0.8;
                	double minDistance = 2*(innerRadius*innerRadius)*(1+Math.cos(etomica.units.Degree.UNIT.toSim(27.0)));
                	if (dr.squared() < minDistance){
                		return 0;
                	}
            	}
        		previousAtom = thisAtom;
        		thisAtom = nextAtom;
        	}
    	}
    	return 1;
    }
    /**
     * Returns log of the ratio of the trial probabilities, ln(Tij/Tji) for the
     * states encountered before (i) and after (j) the most recent call to doTrial(). 
     * Tij is the probability that this move would generate state j from state i, and
     * Tji is the probability that a subsequent call to doTrial would return to state i
     * from state j.
     */
    public double getA() {
    	if (populateList(newSmerList)== 0){
    		return 0;
    	}
    	if (smerList.getAtomCount()!= newSmerList.getAtomCount()){
			return 0;
		}
    	return 1.0;
	}
    
    /**
     * Returns the log of the limiting-distribution probabilities of states, ln(Pj/Pi), 
     * for the states encountered before (i) and after (j) the most recent call to 
     * doTrial.
     */
    public double getB() {
        uNew = 0.0;
        for (int i=0; i<smerList.getAtomCount(); i+=1){
        	energyMeter.setTarget(smerList.getAtom(i));
            uNew += energyMeter.getDataAsScalar();
        }
        return -(uNew - uOld);
    }
    
    public double energyChange() {return uNew - uOld;}
    
    /**
     * Method called by IntegratorMC in the event that the most recent trial is accepted.
     */
    public void acceptNotify() {  /* do nothing */
//    	if (smerList.indexOf(box.getLeafList().getAtom(388))>-1 ||smerList.indexOf(box.getLeafList().getAtom(115))>-1 ){
//        	System.out.println("accepted moving smerList "+smerList);
//        }
    }
    
    /**
     * Method called by IntegratorMC in the event that the most recent trial move is
     * rejected.  This method should cause the system to be restored to the condition
     * before the most recent call to doTrial.
     */
    public void rejectNotify() {
        translationVector.TE(-1);
        for (int i=0; i<smerList.getAtomCount(); i+=1){
        	(smerList.getAtom(i)).getPosition().PE(translationVector);
        	
        }
//        if (smerList.indexOf(box.getLeafList().getAtom(388))>-1 ||smerList.indexOf(box.getLeafList().getAtom(115))>-1 ){
//        	System.out.println("rejected moving smerList "+smerList);
//        	System.out.println("position 388 "+((IAtomPositioned)box.getLeafList().getAtom(388)).getPosition()+"position 115 "+((IAtomPositioned)box.getLeafList().getAtom(115)).getPosition());
//        }
    }
        
    
    public AtomIterator affectedAtoms() {
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
     * @param atomSource The atomSource to set.
     */
    public void setAtomSource(AtomSource source) {
        atomSource = source;
    }
}
