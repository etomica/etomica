package etomica.paracetamol;

/*
 * The agent source for the charges on paracetamol
 */

import etomica.atom.IAtom;
import etomica.atom.IMolecule;
import etomica.atom.AtomAgentManager.AgentSource;
import etomica.potential.EwaldSummation.MyCharge;

public class ChargeAgentSourceParacetamol implements AgentSource {
	
	public ChargeAgentSourceParacetamol(){
		myCharge = new MyCharge[AtomParacetamol.Echarge.length];
		
		for (int i=0; i<AtomParacetamol.Echarge.length; i++){
			myCharge[i] = new MyCharge(AtomParacetamol.Echarge[i]);
			
		}
	}
	
	public Class getAgentClass() {
		
		return MyCharge.class;
	}

	public Object makeAgent(IAtom a) {
		
		if (a instanceof IMolecule){
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
