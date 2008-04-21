package etomica.models.water;

import etomica.api.IAtomLeaf;
import etomica.atom.AtomLeafAgentManager;
import etomica.potential.EwaldSummation.MyCharge;

public class ChargeAgentSourceTIP4PWater implements AtomLeafAgentManager.AgentSource {
	
	public ChargeAgentSourceTIP4PWater(){
		myCharge = new MyCharge[ConformationWaterTIP4P.Echarge.length];
		
		for (int i=0; i<ConformationWaterTIP4P.Echarge.length; i++){
			myCharge[i] = new MyCharge(ConformationWaterTIP4P.Echarge[i]);
			
		}
	}
	
	public Class getAgentClass() {
		
		return MyCharge.class;
	}

	public Object makeAgent(IAtomLeaf a) {
		
		int index = a.getIndex();
		
		return myCharge[index];
	}

	public void releaseAgent(Object agent, IAtomLeaf atom) {
		// Do nothing

	}
	
	protected final MyCharge[] myCharge;
}
