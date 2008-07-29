package etomica.potential;

import etomica.api.IAtomLeaf;
import etomica.api.IAtomSet;
import etomica.api.IMolecule;
import etomica.api.IPotential;
import etomica.api.IVector;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.MoleculeAgentManager;
import etomica.integrator.Integrator;
import etomica.integrator.IntegratorBox;

/**
 * Sums the force and torque on each iterated atom or molecule and adds it to
 * the agent associated with the atom.
 */
public class PotentialCalculationTorqueSum implements PotentialCalculation {
        
    private static final long serialVersionUID = 1L;
    protected AtomLeafAgentManager leafAgentManager;
    protected AtomLeafAgentManager.AgentIterator leafAgentIterator;
    protected MoleculeAgentManager moleculeAgentManager;
    protected MoleculeAgentManager.AgentIterator moleculeAgentIterator;
    
    public void setAgentManager(AtomLeafAgentManager agentManager) {
        leafAgentManager = agentManager;
        leafAgentIterator = leafAgentManager.makeIterator();
    }
    
    public void setMoleculeAgentManager(MoleculeAgentManager newMoleculeAgentManager) {
        moleculeAgentManager = newMoleculeAgentManager;
        moleculeAgentIterator = moleculeAgentManager.makeIterator();
    }
    
    /**
     * Re-zeros the force vectors.
     *
     */
    public void reset(){
        
        leafAgentIterator.reset();
        while(leafAgentIterator.hasNext()){
            Object agent = leafAgentIterator.next();
            if (agent instanceof Integrator.Torquable) {
                ((Integrator.Torquable)agent).torque().E(0);
                ((Integrator.Forcible)agent).force().E(0);
            }
            else if (agent instanceof Integrator.Forcible) {
                ((Integrator.Forcible)agent).force().E(0);
            }
        }
        moleculeAgentIterator.reset();
        while(moleculeAgentIterator.hasNext()){
            Object agent = moleculeAgentIterator.next();
            if (agent instanceof Integrator.Torquable) {
                ((Integrator.Torquable)agent).torque().E(0);
                ((Integrator.Forcible)agent).force().E(0);
            }
            else if (agent instanceof Integrator.Forcible) {
                ((Integrator.Forcible)agent).force().E(0);
            }
        }

    }
    
    /**
	 * Adds forces and torques due to given potential acting on the atoms produced by the iterator.
	 * Implemented for 1-, 2- and N-body potentials.
	 */
	public void doCalculation(IAtomSet atoms, IPotential potential) {
        int nBody = potential.nBody();
        if (atoms.getAtom(0) instanceof IMolecule) {
            if (potential instanceof IPotentialTorque) {
                // IPotentialTorque will give us gradient and torque in one call
                IPotentialTorque potentialSoft = (IPotentialTorque)potential;
                IVector[][] gt = potentialSoft.gradientAndTorque(atoms);
                IVector[] g = gt[0];
                IVector[] t = gt[1];
                switch(nBody) {
                    case 1:
                        ((IntegratorBox.Torquable)moleculeAgentManager.getAgent((IMolecule)atoms.getAtom(0))).torque().PE(t[0]);
                        ((IntegratorBox.Forcible)moleculeAgentManager.getAgent((IMolecule)atoms.getAtom(0))).force().ME(g[0]);
                        break;
                    case 2:
                        ((IntegratorBox.Torquable)moleculeAgentManager.getAgent((IMolecule)atoms.getAtom(0))).torque().PE(t[0]);
                        ((IntegratorBox.Torquable)moleculeAgentManager.getAgent((IMolecule)atoms.getAtom(1))).torque().PE(t[1]);
                        ((IntegratorBox.Forcible)moleculeAgentManager.getAgent((IMolecule)atoms.getAtom(0))).force().ME(g[0]);
                        ((IntegratorBox.Forcible)moleculeAgentManager.getAgent((IMolecule)atoms.getAtom(1))).force().ME(g[1]);
                        break;
                    default:
                        //XXX atoms.count might not equal f.length.  The potential might size its 
                        //array of vectors to be large enough for one AtomSet and then not resize it
                        //back down for another AtomSet with fewer atoms.
                        for (int i=0; i<atoms.getAtomCount(); i++) {
                            ((IntegratorBox.Torquable)moleculeAgentManager.getAgent((IMolecule)atoms.getAtom(i))).torque().PE(t[i]);
                            ((IntegratorBox.Forcible)moleculeAgentManager.getAgent((IMolecule)atoms.getAtom(i))).force().ME(g[i]);
                        }
                }
            }
            else if (potential instanceof PotentialSoft) {
                // we can only get the gradient, but we're probably just dealing with a set of (leaf) Atoms.
                PotentialSoft potentialSoft = (PotentialSoft)potential;
                IVector[] gradient = potentialSoft.gradient(atoms);
                switch(nBody) {
                    case 1:
                        ((IntegratorBox.Forcible)leafAgentManager.getAgent((IAtomLeaf)atoms.getAtom(0))).force().ME(gradient[0]);
                        break;
                    case 2:
                        ((IntegratorBox.Forcible)leafAgentManager.getAgent((IAtomLeaf)atoms.getAtom(0))).force().ME(gradient[0]);
                        ((IntegratorBox.Forcible)leafAgentManager.getAgent((IAtomLeaf)atoms.getAtom(1))).force().ME(gradient[1]);
                        break;
                    default:
                        //XXX atoms.count might not equal f.length.  The potential might size its 
                        //array of vectors to be large enough for one AtomSet and then not resize it
                        //back down for another AtomSet with fewer atoms.
                        for (int i=0; i<atoms.getAtomCount(); i++) {
                            ((IntegratorBox.Forcible)leafAgentManager.getAgent((IAtomLeaf)atoms.getAtom(i))).force().ME(gradient[i]);
                        }
                }
            }
        }
        else {
    	    if (potential instanceof IPotentialTorque) {
    	        // IPotentialTorque will give us gradient and torque in one call
        		IPotentialTorque potentialSoft = (IPotentialTorque)potential;
    			IVector[][] gt = potentialSoft.gradientAndTorque(atoms);
    			IVector[] g = gt[0];
    			IVector[] t = gt[1];
    			switch(nBody) {
    				case 1:
    					((IntegratorBox.Torquable)leafAgentManager.getAgent((IAtomLeaf)atoms.getAtom(0))).torque().PE(t[0]);
                        ((IntegratorBox.Forcible)leafAgentManager.getAgent((IAtomLeaf)atoms.getAtom(0))).force().ME(g[0]);
    					break;
    				case 2:
                        ((IntegratorBox.Torquable)leafAgentManager.getAgent((IAtomLeaf)atoms.getAtom(0))).torque().PE(t[0]);
                        ((IntegratorBox.Torquable)leafAgentManager.getAgent((IAtomLeaf)atoms.getAtom(1))).torque().PE(t[1]);
                        ((IntegratorBox.Forcible)leafAgentManager.getAgent((IAtomLeaf)atoms.getAtom(0))).force().ME(g[0]);
                        ((IntegratorBox.Forcible)leafAgentManager.getAgent((IAtomLeaf)atoms.getAtom(1))).force().ME(g[1]);
    			 		break;
                    default:
                        //XXX atoms.count might not equal f.length.  The potential might size its 
                        //array of vectors to be large enough for one AtomSet and then not resize it
                        //back down for another AtomSet with fewer atoms.
                        for (int i=0; i<atoms.getAtomCount(); i++) {
                            ((IntegratorBox.Torquable)leafAgentManager.getAgent((IAtomLeaf)atoms.getAtom(i))).torque().PE(t[i]);
                            ((IntegratorBox.Forcible)leafAgentManager.getAgent((IAtomLeaf)atoms.getAtom(i))).force().ME(g[i]);
                        }
    			}
    	    }
    	    else if (potential instanceof PotentialSoft) {
                // we can only get the gradient, but we're probably just dealing with a set of (leaf) Atoms.
                PotentialSoft potentialSoft = (PotentialSoft)potential;
                IVector[] gradient = potentialSoft.gradient(atoms);
                switch(nBody) {
                    case 1:
                        ((IntegratorBox.Forcible)leafAgentManager.getAgent((IAtomLeaf)atoms.getAtom(0))).force().ME(gradient[0]);
                        break;
                    case 2:
                        ((IntegratorBox.Forcible)leafAgentManager.getAgent((IAtomLeaf)atoms.getAtom(0))).force().ME(gradient[0]);
                        ((IntegratorBox.Forcible)leafAgentManager.getAgent((IAtomLeaf)atoms.getAtom(1))).force().ME(gradient[1]);
                        break;
                    default:
                        //XXX atoms.count might not equal f.length.  The potential might size its 
                        //array of vectors to be large enough for one AtomSet and then not resize it
                        //back down for another AtomSet with fewer atoms.
                        for (int i=0; i<atoms.getAtomCount(); i++) {
                            ((IntegratorBox.Forcible)leafAgentManager.getAgent((IAtomLeaf)atoms.getAtom(i))).force().ME(gradient[i]);
                        }
                }
    	    }
        }
	}
}
