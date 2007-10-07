package etomica.models.water;

import etomica.atom.IAtom;
import etomica.atom.IAtomGroup;
import etomica.atom.AtomAgentManager.AgentSource;
import etomica.potential.EwaldSummation.MyCharge;

public class ChargeAgentSourceTIP4PWater implements AgentSource {
	
	public ChargeAgentSourceTIP4PWater(){
		myCharge = new MyCharge[AtomWater4P.Echarge.length];
		
		for (int i=0; i<AtomWater4P.Echarge.length; i++){
			myCharge[i] = new MyCharge(AtomWater4P.Echarge[i]);
			
		}
	}
	
	public Class getAgentClass() {
		
		return MyCharge.class;
	}

	public Object makeAgent(IAtom a) {
		
		if (a instanceof IAtomGroup){
			return null;
		}
		
		int index = a.getIndex();
		
		return myCharge[index];
	}

	public void releaseAgent(Object agent, IAtom atom) {
		// Do nothing

	}
	
	protected final MyCharge[] myCharge;
}
