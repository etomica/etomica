package etomica.dimer;

import etomica.atom.AtomSet;
import etomica.atom.IAtomPositioned;
import etomica.atom.iterator.AtomsetIterator;
import etomica.box.Box;
import etomica.integrator.IntegratorBox;
import etomica.potential.IPotential;
import etomica.potential.PotentialCalculationForceSum;
import etomica.potential.PotentialSoft;
import etomica.space.IVector;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.species.ISpecies;

/**
 * Sums the force on each iterated atom and adds it to the integrator agent
 * associated with the atom.  Additionally, this class has the potential
 * calculate the pressureTensor (which can be done efficiently during the
 * gradient calculation).
 */
public class PotentialCalculationForcePressureSumGB extends PotentialCalculationForceSum {
        
    private static final long serialVersionUID = 1L;
    protected final Tensor pressureTensor;
    protected ISpecies fixed;
    protected Box box;
    
    public PotentialCalculationForcePressureSumGB(Space space, Box box) {
        pressureTensor = space.makeTensor();
        this.box = box;
    }
    
    /**
     * Zeros out the pressureTensor.  This method should be called before
     * invoking potentialMaster.calculate so that the pressureTensor is
     * correct at the end of the calculation.
     */
    public void reset() {
        super.reset();
        pressureTensor.E(0);
    }
    
    /**
	 * Adds forces due to given potential acting on the atoms produced by the iterator.
	 * Implemented for only 1- and 2-body potentials.
	 */
	public void doCalculation(AtomsetIterator iterator, IPotential potential) {
		PotentialSoft potentialSoft = (PotentialSoft)potential;
		int nBody = potential.nBody();
		
		IVector forceTop = potential.getSpace().makeVector();
		IVector forceBottom = potential.getSpace().makeVector();
		
		iterator.reset();
		for (AtomSet atoms = iterator.next(); atoms != null; atoms = iterator.next()) {
			IVector[] f = potentialSoft.gradient(atoms, pressureTensor);
			IVector rij = potential.getSpace().makeVector();
			switch(nBody) {
				case 1:
					((IntegratorBox.Forcible)integratorAgentManager.getAgent(atoms.getAtom(0))).force().ME(f[0]);
					break;
				case 2:
                    ((IntegratorBox.Forcible)integratorAgentManager.getAgent(atoms.getAtom(0))).force().ME(f[0]);
                    ((IntegratorBox.Forcible)integratorAgentManager.getAgent(atoms.getAtom(1))).force().ME(f[1]);
			 		break;
                default:
                    //XXX atoms.count might not equal f.length.  The potential might size its 
                    //array of vectors to be large enough for one AtomSet and then not resize it
                    //back down for another AtomSet with fewer atoms.
                    
                    //Find average force in Z-direction and assign to all atoms.
                    for (int i=0; i<atoms.getAtomCount(); i++){
                        rij.E(((IAtomPositioned)atoms.getAtom(i)).getPosition());      
                            if(rij.x(2)>0){
                                forceTop.PE(f[i]);        
                            }
                            else{
                                forceBottom.PE(f[i]);
                            }
                            
                        }
                    } 
			
			        forceTop.TE(2.0/box.atomCount());
			        forceBottom.TE(2.0/box.atomCount());
			        
                    for (int i=0; i<atoms.getAtomCount(); i++){
                        rij.E(((IAtomPositioned)atoms.getAtom(i)).getPosition());
                        
                        if(rij.x(2)>0){
                            f[i].E(forceTop);
                        }
                        else{
                            f[i].E(forceBottom);
                        }
                        ((IntegratorBox.Forcible)integratorAgentManager.getAgent(atoms.getAtom(i))).force().ME(f[i]);
                        
                    }
                    
		}
	}

    /**
     * Returns the virial portion of pressure tensor calculated during the last
     * potential calculation.  In order to be valid, reset() must be called
     * before invoking potentialMaster.calculate.  The given tensor has not
     * been normalized by the system volume.
     */
    public Tensor getPressureTensor() {
        return pressureTensor;
    }
}
