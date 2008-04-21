package etomica.paracetamol;

/*
 * The agent source for the charges on paracetamol
 */

import etomica.api.IAtomLeaf;
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

	public Object makeAgent(IAtomLeaf a) {
		
		int index = a.getIndex();
		
		return myCharge[index];
	}

	public void releaseAgent(Object agent, IAtomLeaf atom) {
		// Do nothing

	}
	
	protected final MyCharge[] myCharge;
}
