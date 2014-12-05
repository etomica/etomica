/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.association;

import etomica.api.IAtom;
import etomica.api.IAtomList;
import etomica.api.IPotentialMaster;
import etomica.api.IRandom;
import etomica.api.ISimulation;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.space.ISpace;
import etomica.virial.simulations.TestLJAssociationMC3D_NPT;

public class MCMoveAtomDimer extends MCMoveAtom {
	protected AssociationManager associationManager;
	
	

	public MCMoveAtomDimer(ISimulation sim, IPotentialMaster potentialMaster,
			ISpace _space) {
		this(potentialMaster, sim.getRandom(), _space, 1.0, 15.0, false);
	}


	public MCMoveAtomDimer(IPotentialMaster potentialMaster, IRandom random,
			ISpace _space, double stepSize, double stepSizeMax,
			boolean fixOverlap) {
		super(potentialMaster, random, _space, stepSize, stepSizeMax,
				fixOverlap);
	}
	public void setAssociationManager(AssociationManager associationManager){
		this.associationManager = associationManager;
		AtomSourceRandomDimer atomSourceRandomDimer = new AtomSourceRandomDimer();
		atomSourceRandomDimer.setAssociationManager(associationManager);
		if (box != null) {
			atomSourceRandomDimer.setBox(box);
		}
		atomSourceRandomDimer.setRandomNumberGenerator(random);
		setAtomSource(atomSourceRandomDimer);
	}
	public double getA(){
		if (associationManager.getAssociatedAtoms(atom).getAtomCount() > 1) {
        	return 0;
        } 
        if (associationManager.getAssociatedAtoms(atom).getAtomCount() == 1){
        	IAtom atomj = associationManager.getAssociatedAtoms(atom).getAtom(0);
        	if(associationManager.getAssociatedAtoms(atomj).getAtomCount() > 1){
        		return 0;
        	} 
        	return 1.0;
        }
    	return 0;
	}

}
