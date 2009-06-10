package etomica.association;

import etomica.api.IAtom;
import etomica.api.IAtomList;
import etomica.api.IPotentialMaster;
import etomica.api.IRandom;
import etomica.api.ISimulation;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.space.ISpace;
import etomica.virial.simulations.TestLJAssociationMC3D_NPT;

public class MCMoveAtomMonomer extends MCMoveAtom {
	protected AssociationManager associationManager;
	
	

	public MCMoveAtomMonomer(ISimulation sim, IPotentialMaster potentialMaster,
			ISpace _space) {
		this(potentialMaster, sim.getRandom(), _space, 1.0, 15.0, false);
	}


	public MCMoveAtomMonomer(IPotentialMaster potentialMaster, IRandom random,
			ISpace _space, double stepSize, double stepSizeMax,
			boolean fixOverlap) {
		super(potentialMaster, random, _space, stepSize, stepSizeMax,
				fixOverlap);
	}
	public void setAssociationManager(AssociationManager associationManager){
		this.associationManager = associationManager;
		AtomSourceRandomMonomer atomSourceRandomMonomer = new AtomSourceRandomMonomer();
		atomSourceRandomMonomer.setAssociationManager(associationManager);
		if (box != null) {
			atomSourceRandomMonomer.setBox(box);
		}
		atomSourceRandomMonomer.setRandomNumberGenerator(random);
		setAtomSource(atomSourceRandomMonomer);
	}

	public double getA(){
		if (associationManager.getAssociatedAtoms(atom).getAtomCount() > 0) {
        	return 0.0;
        } 
		return 1.0;
	}

}
