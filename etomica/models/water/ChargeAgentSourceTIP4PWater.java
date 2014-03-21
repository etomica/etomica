package etomica.models.water;

import etomica.api.IAtom;
import etomica.atom.AtomLeafAgentManager;
import etomica.potential.EwaldSummation.MyCharge;

public class ChargeAgentSourceTIP4PWater implements AtomLeafAgentManager.AgentSource<MyCharge> {
	
	public ChargeAgentSourceTIP4PWater(){
		myCharge = new MyCharge[ConformationWaterTIP4P.Echarge.length];
		
		for (int i=0; i<ConformationWaterTIP4P.Echarge.length; i++){
			myCharge[i] = new MyCharge(ConformationWaterTIP4P.Echarge[i]);
			
		}
	}

	public MyCharge makeAgent(IAtom a) {
		
		int index = a.getIndex();
		
		return myCharge[index];
	}

	public void releaseAgent(MyCharge agent, IAtom atom) {
		// Do nothing

	}
	
	protected final MyCharge[] myCharge;
}
