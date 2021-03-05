/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.association;

import etomica.atom.AtomArrayList;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.potential.PotentialMaster;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.util.random.IRandom;

public class MCMoveAtomMonomer extends MCMoveAtom {
	protected AssociationManager associationManager;
	protected final AtomArrayList smerList;
	protected final Vector dr;
	protected int maxLength = Integer.MAX_VALUE;


	public MCMoveAtomMonomer(IRandom random, PotentialMaster potentialMaster,
							 Space _space) {
		this(potentialMaster, random, _space, 1.0, 15.0, false);
	}


	public MCMoveAtomMonomer(PotentialMaster potentialMaster, IRandom random,
							 Space _space, double stepSize, double stepSizeMax,
							 boolean fixOverlap) {
		super(random, potentialMaster, _space, stepSize, stepSizeMax,
				fixOverlap);
		this.smerList = new AtomArrayList();
		this.dr = _space.makeVector();
	}
	public void setAssociationManager(AssociationManager associationManager){
		this.associationManager = associationManager;
		AtomSourceRandomMonomer atomSourceRandomMonomer = new AtomSourceRandomMonomer();
		atomSourceRandomMonomer.setAssociationManager(associationManager);
		if (box != null) {
			atomSourceRandomMonomer.setBox(box);
		}
		atomSourceRandomMonomer.setRandomNumberGenerator(random);
		setAtomSource(atomSourceRandomMonomer);
	}

    public double getChi(double temperature) {
        if (populateList(smerList)== 0){
    		return 0;
    	}
		if (associationManager.getAssociatedAtoms(atom).size() > 0) {
        	return 0.0;
        } 
		if (smerList.size() > maxLength) {
    		return 0.0;
		}
        return super.getChi(temperature);
    }
	
	public void setMaxLength(int i){
    	maxLength = i;
    }
	
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
    		dr.Ev1Mv2((atom0).getPosition(), (atom1).getPosition());//dr = distance from the atom0 to atom1
        	box.getBoundary().nearestImage(dr);
        	double innerRadius = 0.8;
        	double minDistance = 2*(innerRadius*innerRadius)*(1+Math.cos(etomica.units.Degree.UNIT.toSim(27.0)));
        	if (dr.squared() < minDistance){
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
    		dr.Ev1Mv2((atom0).getPosition(), (atom1).getPosition());//dr = distance from the atom0 to atom1
        	box.getBoundary().nearestImage(dr);
        	double innerRadius = 0.8;
        	double minDistance = 2*(innerRadius*innerRadius)*(1+Math.cos(etomica.units.Degree.UNIT.toSim(27.0)));
        	if (dr.squared() < minDistance){
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
        		dr.Ev1Mv2((atom0).getPosition(), (atom1).getPosition());//dr = distance from the atom0 to atom1
            	box.getBoundary().nearestImage(dr);
            	double innerRadius = 0.8;
            	double minDistance = 2*(innerRadius*innerRadius)*(1+Math.cos(etomica.units.Degree.UNIT.toSim(27.0)));
            	if (dr.squared() < minDistance){
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
}
