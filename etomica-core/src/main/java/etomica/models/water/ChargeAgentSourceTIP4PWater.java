/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.water;

import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.atom.AtomLeafAgentManager;
import etomica.potential.EwaldSummation.MyCharge;

public class ChargeAgentSourceTIP4PWater implements AtomLeafAgentManager.AgentSource<MyCharge> {
	
	public ChargeAgentSourceTIP4PWater(){
		myCharge = new MyCharge[ConformationWaterTIP4P.Echarge.length];
		
		for (int i=0; i<ConformationWaterTIP4P.Echarge.length; i++){
			myCharge[i] = new MyCharge(ConformationWaterTIP4P.Echarge[i]);
			
		}
	}

	public MyCharge makeAgent(IAtom a, Box agentBox) {
		
		int index = a.getIndex();
		
		return myCharge[index];
	}

	public void releaseAgent(MyCharge agent, IAtom atom, Box agentBox) {
		// Do nothing

	}
	
	protected final MyCharge[] myCharge;
}
