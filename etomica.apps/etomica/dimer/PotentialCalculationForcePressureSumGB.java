package etomica.dimer;

import etomica.atom.AtomSet;
import etomica.atom.IAtomPositioned;
import etomica.atom.IMolecule;
import etomica.atom.iterator.AtomsetIterator;
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
    
    public PotentialCalculationForcePressureSumGB(Space space, ISpecies fixedSpecies) {
        pressureTensor = space.makeTensor();
        this.fixed = fixedSpecies;
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
		double fixZsum = 0.0;
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
                    
                    //Find average force in Z-direction and assign to all fixed atoms.
                    for (int i=0; i<atoms.getAtomCount(); i++){
                        if(atoms.getAtom(i).getType().getSpecies()==fixed){
                            fixZsum += f[i].x(2) / fixed.getNumLeafAtoms();
                        }
                    }
                    for (int i=0; i<atoms.getAtomCount(); i++){
                        if(atoms.getAtom(i).getType().getSpecies()==fixed){
                            rij.E(((IAtomPositioned)((IMolecule)atoms.getAtom(i)).getChildList().getAtom(0)).getPosition());
                            if(rij.x(2) < 0){
                                f[i].setX(2,-fixZsum);
                            }
                            else{
                                f[i].setX(2,fixZsum);
                            }
                        }
                    }
                    
                    for (int i=0; i<atoms.getAtomCount(); i++) {
                        ((IntegratorBox.Forcible)integratorAgentManager.getAgent(atoms.getAtom(i))).force().ME(f[i]);
                    }
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
