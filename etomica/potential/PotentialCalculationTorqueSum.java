package etomica.potential;

import etomica.api.IVector;
import etomica.atom.AtomAgentManager;
import etomica.atom.AtomSet;
import etomica.atom.iterator.AtomsetIterator;
import etomica.integrator.IntegratorBox;

/**
 * Sums the force and torque on each iterated atom or molecule and adds it to
 * the agent associated with the atom.
 */
public class PotentialCalculationTorqueSum extends PotentialCalculation {
        
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
            Object agent = agentIterator.next();
            if (agent instanceof IntegratorBox.Torquable) {
                ((IntegratorBox.Torquable)agent).torque().E(0);
                ((IntegratorBox.Forcible)agent).force().E(0);
            }
            else if (agent instanceof IntegratorBox.Forcible) {
                ((IntegratorBox.Forcible)agent).force().E(0);
            }
        }
    }
    
    /**
	 * Adds forces and torques due to given potential acting on the atoms produced by the iterator.
	 * Implemented for 1-, 2- and N-body potentials.
	 */
	public void doCalculation(AtomsetIterator iterator, IPotential potential) {
	    if (potential instanceof IPotentialTorque) {
	        // IPotentialTorque will give us gradient and torque in one call
    		IPotentialTorque potentialSoft = (IPotentialTorque)potential;
    		int nBody = potential.nBody();
    		iterator.reset();
            for (AtomSet atoms = iterator.next(); atoms != null; atoms = iterator.next()) {
    			IVector[][] gt = potentialSoft.gradientAndTorque(atoms);
    			IVector[] g = gt[0];
    			IVector[] t = gt[1];
    			switch(nBody) {
    				case 1:
    					((IntegratorBox.Torquable)integratorAgentManager.getAgent(atoms.getAtom(0))).torque().PE(t[0]);
                        ((IntegratorBox.Forcible)integratorAgentManager.getAgent(atoms.getAtom(0))).force().ME(g[0]);
    					break;
    				case 2:
                        ((IntegratorBox.Torquable)integratorAgentManager.getAgent(atoms.getAtom(0))).torque().PE(t[0]);
                        ((IntegratorBox.Torquable)integratorAgentManager.getAgent(atoms.getAtom(1))).torque().PE(t[1]);
                        ((IntegratorBox.Forcible)integratorAgentManager.getAgent(atoms.getAtom(0))).force().ME(g[0]);
                        ((IntegratorBox.Forcible)integratorAgentManager.getAgent(atoms.getAtom(1))).force().ME(g[1]);
    			 		break;
                    default:
                        //XXX atoms.count might not equal f.length.  The potential might size its 
                        //array of vectors to be large enough for one AtomSet and then not resize it
                        //back down for another AtomSet with fewer atoms.
                        for (int i=0; i<atoms.getAtomCount(); i++) {
                            ((IntegratorBox.Torquable)integratorAgentManager.getAgent(atoms.getAtom(i))).torque().PE(t[i]);
                            ((IntegratorBox.Forcible)integratorAgentManager.getAgent(atoms.getAtom(i))).force().ME(g[i]);
                        }
    			}
    		}
	    }
	    else if (potential instanceof PotentialSoft) {
            // we can only get the gradient, but we're probably just dealing with a set of (leaf) Atoms.
            PotentialSoft potentialSoft = (PotentialSoft)potential;
            int nBody = potential.nBody();
            iterator.reset();
            for (AtomSet atoms = iterator.next(); atoms != null; atoms = iterator.next()) {
                IVector[] gradient = potentialSoft.gradient(atoms);
                switch(nBody) {
                    case 1:
                        ((IntegratorBox.Forcible)integratorAgentManager.getAgent(atoms.getAtom(0))).force().ME(gradient[0]);
                        break;
                    case 2:
                        ((IntegratorBox.Forcible)integratorAgentManager.getAgent(atoms.getAtom(0))).force().ME(gradient[0]);
                        ((IntegratorBox.Forcible)integratorAgentManager.getAgent(atoms.getAtom(1))).force().ME(gradient[1]);
                        break;
                    default:
                        //XXX atoms.count might not equal f.length.  The potential might size its 
                        //array of vectors to be large enough for one AtomSet and then not resize it
                        //back down for another AtomSet with fewer atoms.
                        for (int i=0; i<atoms.getAtomCount(); i++) {
                            ((IntegratorBox.Forcible)integratorAgentManager.getAgent(atoms.getAtom(i))).force().ME(gradient[i]);
                        }
                }
            }
	    }
	}
}
