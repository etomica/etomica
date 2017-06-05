/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.paracetamol;

/*
 * The agent source for the charges on paracetamol
 */

import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.atom.AtomLeafAgentManager.AgentSource;
import etomica.potential.EwaldSummation.MyCharge;

public class ChargeAgentSourceParacetamol implements AgentSource {
	
	public ChargeAgentSourceParacetamol(){
		myCharge = new MyCharge[SpeciesParacetamol.Echarge.length];
		
		for (int i=0; i<SpeciesParacetamol.Echarge.length; i++){
			myCharge[i] = new MyCharge(SpeciesParacetamol.Echarge[i]);
			
		}
	}
	
	public Class getAgentClass() {
		
		return MyCharge.class;
	}

	public Object makeAgent(IAtom a, Box agentBox) {
		
		int index = a.getIndex();
		
		return myCharge[index];
	}

	public void releaseAgent(Object agent, IAtom atom, Box agentBox) {
		// Do nothing

	}
	
	protected final MyCharge[] myCharge;
}
