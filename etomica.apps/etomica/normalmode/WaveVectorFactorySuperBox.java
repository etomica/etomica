package etomica.normalmode;

import etomica.api.IBox;
import etomica.api.IVectorMutable;
import etomica.lattice.crystal.Primitive;
import etomica.space.ISpace;

public class WaveVectorFactorySuperBox extends WaveVectorFactorySimple {

	/**
	 * Specifically for FCC 32 system size
	 * 
	 * Instead of making wavevectors for 216 cells; this class will
	 * 	reduce the size of box's boundary by one-third in order 
	 * 	to make wave vectors of 8 cells
	 * 
	 * and then the size of box's boundary is restored to the initial
	 * 	condition. 
	 * 
	 * @author Tai Boon Tan
	 */
	private static final long serialVersionUID = 1L;

	public WaveVectorFactorySuperBox(Primitive primitive, ISpace _space) {
		super(primitive, _space);
		// TODO Auto-generated constructor stub
	}
	
	public void makeWaveVectors(IBox box){
		
		IVectorMutable boxDimension = space.makeVector();
		boxDimension.E(box.getBoundary().getDimensions());
		double fraction = (double)(1.0/3);
		boxDimension.TE(fraction);
		box.getBoundary().setDimensions(boxDimension);
		
		super.makeWaveVectors(box);
		
		boxDimension.TE(3.0);
		box.getBoundary().setDimensions(boxDimension);
		
	}
	
}
