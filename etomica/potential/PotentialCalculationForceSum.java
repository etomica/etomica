package etomica.potential;

import etomica.atom.AtomAgentManager;
import etomica.atom.AtomPair;
import etomica.atom.AtomSet;
import etomica.atom.IAtom;
import etomica.atom.iterator.AtomsetIterator;
import etomica.integrator.IntegratorPhase;
import etomica.space.IVector;

/**
 * Sums the force on each iterated atom and adds it to the integrator agent
 * associated with the atom.
 */
public class PotentialCalculationForceSum extends PotentialCalculation {
        
    private static final long serialVersionUID = 1L;
    protected AtomAgentManager integratorAgentManager;
    protected AtomAgentManager.AgentIterator agentIterator;
    
    public void setAgentManager(AtomAgentManager agentManager) {
        integratorAgentManager = agentManager;
        agentIterator = integratorAgentManager.makeIterator();
    }
    /**
     * Re-zeros the force vectors.
     *
     */
    public void reset(){
        
        agentIterator.reset();
        while(agentIterator.hasNext()){
            ((IntegratorPhase.Forcible)agentIterator.next()).force().E(0);
        }
    }
    
    /**
	 * Adds forces due to given potential acting on the atoms produced by the iterator.
	 * Implemented for only 1- and 2-body potentials.
	 */
	public void doCalculation(AtomsetIterator iterator, Potential potential) {
		PotentialSoft potentialSoft = (PotentialSoft)potential;
		int nBody = potential.nBody();
		iterator.reset();
		while(iterator.hasNext()) {
			AtomSet atoms = iterator.next();
			IVector[] f = potentialSoft.gradient(atoms);
			switch(nBody) {
				case 1:
					((IntegratorPhase.Forcible)integratorAgentManager.getAgent((IAtom)atoms)).force().ME(f[0]);
					break;
				case 2:
                    ((IntegratorPhase.Forcible)integratorAgentManager.getAgent(((AtomPair)atoms).atom0)).force().ME(f[0]);
                    ((IntegratorPhase.Forcible)integratorAgentManager.getAgent(((AtomPair)atoms).atom1)).force().ME(f[1]);
			 		break;
                default:
                    //XXX atoms.count might not equal f.length.  The potential might size its 
                    //array of vectors to be large enough for one AtomSet and then not resize it
                    //back down for another AtomSet with fewer atoms.
                    for (int i=0; i<atoms.count(); i++) {
                        ((IntegratorPhase.Forcible)integratorAgentManager.getAgent(atoms.getAtom(i))).force().ME(f[i]);
                    }
			}
		}
	}
}
