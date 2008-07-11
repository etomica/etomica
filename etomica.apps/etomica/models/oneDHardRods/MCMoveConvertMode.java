package etomica.models.oneDHardRods;

import etomica.api.IBox;
import etomica.api.IPotentialMaster;
import etomica.atom.iterator.AtomIterator;
import etomica.integrator.mcmove.MCMove;
import etomica.normalmode.CoordinateDefinition;

public class MCMoveConvertMode extends MCMove {

	CoordinateDefinition cd;
	double temperature;
	IBox box;
	
	
	public MCMoveConvertMode(IPotentialMaster potentialMaster) {
		super(potentialMaster);
	}

	@Override
	public void acceptNotify() {
		// TODO Auto-generated method stub
		
	}

	@Override
	public AtomIterator affectedAtoms(IBox box) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public boolean doTrial() {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public double energyChange(IBox box) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public double getA() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public double getB() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public void rejectNotify() {
		// TODO Auto-generated method stub
		
	}

	void setCoordinateDefinition(CoordinateDefinition cd){
		this.cd = cd;
	}
	
	void setBox(IBox box){
		this.box = box;
	}
}
