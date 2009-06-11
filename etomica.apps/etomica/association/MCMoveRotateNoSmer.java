package etomica.association;

import etomica.api.IPotentialMaster;
import etomica.api.IRandom;
import etomica.integrator.mcmove.MCMoveRotate;
import etomica.space.ISpace;

public class MCMoveRotateNoSmer extends MCMoveRotate {

	public MCMoveRotateNoSmer(IPotentialMaster potentialMaster, IRandom random,
			ISpace _space) {
		super(potentialMaster, random, _space);
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
