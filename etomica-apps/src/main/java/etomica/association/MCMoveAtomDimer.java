/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.association;

import etomica.atom.IAtom;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.util.random.IRandom;

public class MCMoveAtomDimer extends MCMoveAtom {
	protected AssociationManager associationManager;
	
	

	public MCMoveAtomDimer(Simulation sim, PotentialMaster potentialMaster,
                           Space _space) {
		this(potentialMaster, sim.getRandom(), _space, 1.0, 15.0, false);
	}


	public MCMoveAtomDimer(PotentialMaster potentialMaster, IRandom random,
                           Space _space, double stepSize, double stepSizeMax,
                           boolean fixOverlap) {
		super(random, potentialMaster, _space, stepSize, stepSizeMax,
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

	public double getChi(double temperature) {
		if (associationManager.getAssociatedAtoms(atom).size() > 1) {
        	return 0;
        } 
        if (associationManager.getAssociatedAtoms(atom).size() == 1){
        	IAtom atomj = associationManager.getAssociatedAtoms(atom).get(0);
        	if(associationManager.getAssociatedAtoms(atomj).size() > 1){
        		return 0;
        	}
			return super.getChi(temperature);
		}
    	return 0;
	}

}
