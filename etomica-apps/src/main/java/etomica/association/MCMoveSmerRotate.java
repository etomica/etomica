/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.association;

import etomica.atom.*;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorArrayListSimple;
import etomica.box.Box;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.mcmove.MCMoveBoxStep;
import etomica.nbr.cell.Api1ACell;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.potential.IPotentialAtomic;
import etomica.simulation.Simulation;
import etomica.space.RotationTensor;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.util.random.IRandom;

/**
 * Standard Monte Carlo atom-displacement trial move.
 *
 * @author David Kofke
 */
public class MCMoveSmerRotate extends MCMoveBoxStep {
    
    private static final long serialVersionUID = 2L;
    protected final AtomIteratorArrayListSimple affectedAtomIterator;
    protected final MeterPotentialEnergy energyMeter;
    protected IAtom atom;
    protected double uOld;
    protected double uNew = Double.NaN;
    protected AtomSource atomSource;
    protected boolean fixOverlap;
    protected final IRandom random;
    protected Space space;
    protected final PotentialMasterCell potentialMaster;
    protected final Api1ACell neighborIterator;
    protected final Vector dr, dr2;
    protected final IPotentialAtomic trimerPotential;
    protected final AtomArrayList smerList;
    protected final AtomArrayList newSmerList;
    protected transient RotationTensor rotationTensor;
    protected final Vector r0;
    protected AssociationManager associationManager;

    public MCMoveSmerRotate(Simulation sim, PotentialMasterCell potentialMaster, Space _space, IPotentialAtomic trimerPotential) {
        this(potentialMaster, sim.getRandom(), _space, 1.0, 15.0, false, trimerPotential);
    }
    
    public MCMoveSmerRotate(PotentialMasterCell potentialMaster, IRandom random,
                            Space _space, double stepSize, double stepSizeMax,
                            boolean fixOverlap, IPotentialAtomic trimerPotential) {
        super(potentialMaster);
        this.smerList = new AtomArrayList();
        this.newSmerList = new AtomArrayList();
        this.affectedAtomIterator = new AtomIteratorArrayListSimple(smerList);
        this.potentialMaster = potentialMaster;
        this.neighborIterator = new Api1ACell(3,1.0,potentialMaster.getCellAgentManager());
        this.dr = _space.makeVector();
        this.dr2 = _space.makeVector();
        rotationTensor = _space.makeRotationTensor();
        this.r0= _space.makeVector();
        this.random = random;
        this.space = _space;
        this.trimerPotential = trimerPotential;
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
        populateList(smerList);
        uOld = 0.0;
        for (int i = 0; i<smerList.size(); i+=1){
        	energyMeter.setTarget(smerList.get(i));
            uOld += energyMeter.getDataAsScalar();
        }  
        if(uOld > 1e8 && !fixOverlap) {
            throw new RuntimeException("atom "+atom+" in box "+box+" has an overlap");
        }
        double dTheta = (2*random.nextDouble() - 1.0)*stepSize;
        rotationTensor.setAxial(r0.getD() == 3 ? random.nextInt(3) : 2,dTheta);
        r0.E(0.0);
        for (int i = 0; i<smerList.size(); i+=1){
        	dr.Ev1Mv2(smerList.get(i).getPosition(), r0);//dr = distance from the atom to the center of the mass
        	box.getBoundary().nearestImage(dr);
        	r0.PEa1Tv1(1.0/(i+1), dr);//
        }
        for (int i = 0; i<smerList.size(); i+=1){
        	doTransform((IAtomOriented)smerList.get(i));
        }
    
//        if (atom.getParentGroup().getIndex() == 371 || atom.getParentGroup().getIndex() == 224 ||
//        		((IAtom)atom1).getParentGroup().getIndex() == 371 || ((IAtom)atom1).getParentGroup().getIndex() == 224){
//        	System.out.println("MCMoveDimerRotate "+atom);
//        }
        
        return true;
    }//end of doTrial
    protected int populateList(AtomArrayList mySmerList){
    	mySmerList.clear();
    	mySmerList.add(atom);
    	IAtomList bondList = associationManager.getAssociatedAtoms(atom);
    	if (bondList.size() > 2){
    		return 0;
    	}
    	if (bondList.size() == 2){
    		IAtom atom0 = bondList.get(0);
    		IAtom atom1 = bondList.get(1);
    		dr2.Ev1Mv2((atom0).getPosition(), (atom1).getPosition());
        	box.getBoundary().nearestImage(dr2);
        	double innerRadius = 0.8;
        	double minDistance = 2*(innerRadius*innerRadius)*(1+Math.cos(etomica.units.Degree.UNIT.toSim(27.0)));
        	if (dr2.squared() < minDistance){
        		return 0;
        	}
    	}
    	if (bondList.size() == 0){
    		return 1;
    	}
    	IAtom thisAtom = bondList.get(0);
    	mySmerList.add(thisAtom);
    	IAtomList bondList1 = associationManager.getAssociatedAtoms(thisAtom);
    	if (bondList1.size() > 2){
    		return 0;
    	}
    	if (bondList1.size() == 2){
    		IAtom atom0 = bondList1.get(0);
    		IAtom atom1 = bondList1.get(1);
    		dr2.Ev1Mv2((atom0).getPosition(), (atom1).getPosition());//dr = distance from the atom0 to atom1
        	box.getBoundary().nearestImage(dr2);
        	double innerRadius = 0.8;
        	double minDistance = 2*(innerRadius*innerRadius)*(1+Math.cos(etomica.units.Degree.UNIT.toSim(27.0)));
        	if (dr2.squared() < minDistance){
        		return 0;
        	}
    	}
    	IAtom previousAtom = atom;
    	while (bondList1.size() > 1){
    		IAtom nextAtom = bondList1.get(0);
    		if (nextAtom == previousAtom){
    			nextAtom = bondList1.get(1);
    		} 
    		if (nextAtom == atom){
    			return 1;
    		}
    		mySmerList.add(nextAtom);
    		bondList1 = associationManager.getAssociatedAtoms(nextAtom);
    		if (bondList1.size() > 2){
        		return 0;
        	}
    		if (bondList1.size() == 2){
        		IAtom atom0 = bondList1.get(0);
        		IAtom atom1 = bondList1.get(1);
        		dr2.Ev1Mv2((atom0).getPosition(), (atom1).getPosition());
            	box.getBoundary().nearestImage(dr2);
            	double innerRadius = 0.8;
            	double minDistance = 2*(innerRadius*innerRadius)*(1+Math.cos(etomica.units.Degree.UNIT.toSim(27.0)));
            	if (dr2.squared() < minDistance){
            		return 0;
            	}
        	}
    		previousAtom = thisAtom;
    		thisAtom = nextAtom;
    	}
    	if (bondList.size()>1){
    		thisAtom = bondList.get(1);
        	mySmerList.add(thisAtom);
        	bondList1 = associationManager.getAssociatedAtoms(thisAtom);
        	if (bondList1.size() > 2){
        		return 0;
        	}
        	if (bondList1.size() == 2){
        		IAtom atom0 = bondList1.get(0);
        		IAtom atom1 = bondList1.get(1);
        		dr2.Ev1Mv2((atom0).getPosition(), (atom1).getPosition());
            	box.getBoundary().nearestImage(dr2);
            	double innerRadius = 0.8;
            	double minDistance = 2*(innerRadius*innerRadius)*(1+Math.cos(etomica.units.Degree.UNIT.toSim(27.0)));
            	if (dr2.squared() < minDistance){
            		return 0;
            	}
        	}
        	previousAtom = atom;
        	while (bondList1.size() > 1){
        		IAtom nextAtom = bondList1.get(0);
        		if (nextAtom == previousAtom){
        			nextAtom = bondList1.get(1);
        		} 
        		mySmerList.add(nextAtom);
        		bondList1 = associationManager.getAssociatedAtoms(nextAtom);
        		if (bondList1.size() > 2){
            		return 0;
            	}
        		if (bondList1.size() == 2){
            		IAtom atom0 = bondList1.get(0);
            		IAtom atom1 = bondList1.get(1);
            		dr2.Ev1Mv2((atom0).getPosition(), (atom1).getPosition());
                	box.getBoundary().nearestImage(dr2);
                	double innerRadius = 0.8;
                	double minDistance = 2*(innerRadius*innerRadius)*(1+Math.cos(etomica.units.Degree.UNIT.toSim(27.0)));
                	if (dr2.squared() < minDistance){
                		return 0;
                	}
            	}
        		previousAtom = thisAtom;
        		thisAtom = nextAtom;
        	}
    	}
    	return 1;
    }
    protected void doTransform(IAtomOriented a) {
            Vector r = a.getPosition();
            r.ME(r0);
            box.getBoundary().nearestImage(r);
            rotationTensor.transform(r);
            r.PE(r0);
            dr.E(a.getOrientation().getDirection());
            rotationTensor.transform(dr);
            a.getOrientation().setDirection(dr);
            
        }

    public double getChi(double temperature) {
        if (populateList(newSmerList) == 0) {
            return 0;
        }
        if (smerList.size() != newSmerList.size()) {
            return 0;
		}

    	uNew = 0.0;
        for (int i = 0; i<smerList.size(); i+=1){
        	energyMeter.setTarget(smerList.get(i));
            uNew += energyMeter.getDataAsScalar();
        }
        return Math.exp(-(uNew - uOld) / temperature);
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
        for (int i = 0; i<smerList.size(); i+=1){
        	doTransform((IAtomOriented)smerList.get(i));
        }
    }
        
    
    public AtomIterator affectedAtoms() {;
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
