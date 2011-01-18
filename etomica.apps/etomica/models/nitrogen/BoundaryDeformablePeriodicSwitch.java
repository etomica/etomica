package etomica.models.nitrogen;

import etomica.api.IVector;
import etomica.api.IVectorMutable;
import etomica.space.BoundaryDeformablePeriodic;
import etomica.space.ISpace;

public class BoundaryDeformablePeriodicSwitch extends BoundaryDeformablePeriodic{
	
	public BoundaryDeformablePeriodicSwitch(ISpace _space, IVector[] vex) {
		super(_space, vex);
		temp = space.makeVector();
		
	}
	public void nearestImage(IVectorMutable dr) {
		if(doPBC){
			super.nearestImage(dr);
		} 
	}

	public IVector centralImage(IVector r) {
		
		if(doPBC){
			return super.centralImage(r);
		}
		temp.E(0.0);
		return temp;
	}

	public boolean getPeriodicity(int d) {
		
		return doPBC;
	}

	public boolean isDoPBC() {
		return doPBC;
	}

	public void setDoPBC(boolean doPBC) {
		this.doPBC = doPBC;
	}

	private static final long serialVersionUID = 1L;
	protected IVectorMutable temp;
	protected boolean doPBC = true;
}
