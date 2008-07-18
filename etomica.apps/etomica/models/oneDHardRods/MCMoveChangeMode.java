package etomica.models.oneDHardRods;

import etomica.api.IPotentialMaster;
import etomica.api.IRandom;
import etomica.api.IVector;
import etomica.normalmode.MCMoveHarmonicStep;

public class MCMoveChangeMode extends MCMoveHarmonicStep{

	IVector[] wavevectors;
	
	public MCMoveChangeMode(IPotentialMaster potentialMaster, IRandom random) {
		super(potentialMaster, random);
		 
	
	
	}

	
	public boolean doTrial() {
		//Pick which wavevector will be changed
		
		
		super.doTrial();
		
		return true;
	}


}
