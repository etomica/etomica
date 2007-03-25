package etomica.potential;

import etomica.atom.Atom;
import etomica.atom.AtomPair;
import etomica.atom.AtomSet;
import etomica.atom.iterator.AtomsetIterator;
import etomica.integrator.IntegratorPhase;
import etomica.space.IVector;
import etomica.space.Space;
import etomica.space.Tensor;

/**
 * Sums the force on each iterated atom and adds it to the integrator agent
 * associated with the atom.  Additionally, this class has the potential
 * calculate the pressureTensor (which can be done efficiently during the
 * gradient calculation).
 */
public class PotentialCalculationForcePressureSum extends PotentialCalculationForceSum {
        
    private static final long serialVersionUID = 1L;
    protected final Tensor pressureTensor;
    
    public PotentialCalculationForcePressureSum(Space space) {
        pressureTensor = space.makeTensor();
    }
    
    /**
     * Zeros out the pressureTensor.  This method should be called before
     * invoking potentialMaster.calculate so that the pressureTensor is
     * correct at the end of the calculation.
     */
    public void reset() {
        pressureTensor.E(0);
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
			IVector[] f = potentialSoft.gradient(atoms, pressureTensor);
			switch(nBody) {
				case 1:
					((IntegratorPhase.Forcible)integratorAgentManager.getAgent((Atom)atoms)).force().ME(f[0]);
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

    /**
     * Returns the pressure tensor calculated during the last potential
     * calculation.  In order to be valid, reset() must be called before
     * invoking potentialMaster.calculate.
     */
    public Tensor getPressureTensor() {
        return pressureTensor;
    }
}
