package etomica.association;

import etomica.api.IPotentialMaster;
import etomica.api.IRandom;
import etomica.api.ISimulation;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.space.ISpace;

public class MCMoveAtomNoSmer extends MCMoveAtom {

	public MCMoveAtomNoSmer(ISimulation sim, IPotentialMaster potentialMaster,
			ISpace _space) {
		super(sim, potentialMaster, _space);
	}
	public void setAssociationManager(AssociationManager associationManager){
		this.associationManager = associationManager;
	}
	public double getA() {
        if (associationManager.getAssociatedAtoms(atom).getAtomCount() > 1) {
        	return 0;
        } 
        return 1;
	}
	protected AssociationManager associationManager;
}
