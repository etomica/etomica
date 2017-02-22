/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.dimer;

import etomica.api.IAtomList;
import etomica.api.IBox;
import etomica.api.IPotential;
import etomica.api.ISpecies;
import etomica.api.IVector;
import etomica.api.IVectorMutable;
import etomica.integrator.IntegratorBox;
import etomica.potential.PotentialCalculationForceSum;
import etomica.potential.PotentialSoft;
import etomica.space.ISpace;
import etomica.space.Tensor;


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
    protected IBox box;
    private final ISpace space;
    
    public PotentialCalculationForcePressureSumGB(ISpace _space, IBox box) {
    	this.space = _space;
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
	public void doCalculation(IAtomList atoms, IPotential potential) {
		PotentialSoft potentialSoft = (PotentialSoft)potential;
		int nBody = potential.nBody();
		
		IVectorMutable forceTop = space.makeVector();
		IVectorMutable forceBottom = space.makeVector();
		
		IVector[] f = potentialSoft.gradient(atoms, pressureTensor);
		IVectorMutable rij = space.makeVector();
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
                    rij.E(atoms.getAtom(i).getPosition());      
                        if(rij.getX(2)>0){
                            forceTop.PE(f[i]);        
                        }
                        else{
                            forceBottom.PE(f[i]);
                        }
                }
                //Averages force over all atoms, and subtracts amount from each atom's force.
		        forceTop.TE(2.0/box.getLeafList().getAtomCount());
		        forceBottom.TE(2.0/box.getLeafList().getAtomCount());
		        for (int i=0; i<box.getLeafList().getAtomCount(); i++){
                    rij.E(box.getLeafList().getAtom(i).getPosition());
                    
                    if(rij.getX(2)>0){
                        
                        ((IntegratorBox.Forcible)integratorAgentManager.getAgent(box.getLeafList().getAtom(i))).force().ME(forceTop);
                    }
                    else{
                        
                        ((IntegratorBox.Forcible)integratorAgentManager.getAgent(box.getLeafList().getAtom(i))).force().ME(forceBottom);
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
