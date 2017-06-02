/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.association;

import etomica.api.IAtom;
import etomica.api.IAtomList;
import etomica.potential.PotentialMaster;
import etomica.api.IRandom;
import etomica.space.Vector;
import etomica.atom.AtomArrayList;
import etomica.atom.IAtomOriented;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.space.IOrientation;
import etomica.space.Space;

/**
 * Performs a rotation of an atom (not a molecule) that has an orientation coordinate.
 */
public class MCMoveRotateAssociated extends MCMoveAtom {
    
    private static final long serialVersionUID = 2L;
    private final IOrientation oldOrientation;
    protected AssociationManager associationManager;
    protected final AtomArrayList smerList;
    protected final Vector dr;
	protected int maxLength = Integer.MAX_VALUE;

    private transient IOrientation iOrientation;

    public MCMoveRotateAssociated(PotentialMaster potentialMaster, IRandom random,
                                  Space _space) {
        super(potentialMaster, random, _space, Math.PI/2, Math.PI, false);
        oldOrientation = _space.makeOrientation();
        this.smerList = new AtomArrayList();
        this.dr = _space.makeVector();
    }
    
    public void setAssociationManager(AssociationManager associationManager){
		this.associationManager = associationManager;
    }
    
    public boolean doTrial() {//rotate any atom
        if(box.getMoleculeList().getMoleculeCount()==0) {return false;}
        atom = atomSource.getAtom();

        energyMeter.setTarget(atom);
        uOld = energyMeter.getDataAsScalar();
        iOrientation = ((IAtomOriented)atom).getOrientation(); 
        oldOrientation.E(iOrientation);  //save old orientation
        iOrientation.randomRotation(random, stepSize);
        
        return true;
    }
    
    public double getA(){
    	if (populateList(smerList) == 0){
    		return 0.0;
    	}
    	if (smerList.getAtomCount() > maxLength) {
    		return 0.0;
		}
		return 1.0;
	}
    
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
    		dr.Ev1Mv2((atom0).getPosition(), (atom1).getPosition());
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
    		dr.Ev1Mv2((atom0).getPosition(), (atom1).getPosition());
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
        		dr.Ev1Mv2((atom0).getPosition(), (atom1).getPosition());
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
        		dr.Ev1Mv2((atom0).getPosition(), (atom1).getPosition());
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
            		dr.Ev1Mv2((atom0).getPosition(), (atom1).getPosition());
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
    
    public void setMaxLength(int i){
    	maxLength = i;
    }
    
    public void rejectNotify() {
        iOrientation.E(oldOrientation);
    }
}
