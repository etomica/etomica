/*
 * Created on Jul 7, 2003
 */
package etomica;


import etomica.units.Dimension;
import etomica.utility.*;


/**
 * @author Mitchell Tai
 *
 * Meter which monitors the component k values of a function as a result of
 * a Fourier Transform. 
 */
public class MeterFourierTransform extends MeterFunction {
	AtomIterator iterator; 
	int nMolecules;
	double[] data;
	public boolean REAL=true;		// toggle to meter the real/imaginary transform
	public boolean TRANSFORM=true;	// toggle to transform
	FastFourierTransform fourier;
	

	public MeterFourierTransform() {
		this(Simulation.instance,true);
	}
	public MeterFourierTransform(SimulationElement parent,boolean real) {
		super((Simulation) parent);
		this.REAL=real;
		if (real){setLabel("Real");}
		else {setLabel("Imaginary");}
		setXLabel("Wave Vector (k)");
		fourier = new FastFourierTransform();
	}
	public void setPhase(Phase p) {
		super.setPhase(p);
		iterator = p.makeAtomIterator();
		nPoints = p.atomCount();
		resizeArrays();	
		nMolecules = nPoints;
	}

	public Dimension getXDimension() {
		return Dimension.NULL;
	}
	/**
	 * returns the array of real or imaginary numbers 
	 */
	public double[] getData() {
		iterator.reset();
		Atom currentAtom;
		data = new double[nMolecules];
		for (int i=0;i<nMolecules;i++) {
			currentAtom = iterator.next();
			// calculate the displacement from their home position.
			data[i] = currentAtom.coord.position(0) - (i*Default.BOX_SIZE/nMolecules);
			
			// accomodates for excessive distances due to wrap around
			if (data[i]>Default.BOX_SIZE/2){data[i]-=Default.BOX_SIZE;}
			if (data[i]<-Default.BOX_SIZE/2){data[i]+=Default.BOX_SIZE;}
		}
		fourier.setData(data);
		if(TRANSFORM)fourier.transform();
		x=fourier.getIndex();
		if (REAL){return fourier.getReal();}
		else {return fourier.getImaginary();}
	}

	public Dimension getDimension() {
		return Dimension.NULL;
	}

	public boolean isTRANSFORM() {
		return TRANSFORM;
	}

	public void setTRANSFORM(boolean transform) {
		TRANSFORM = transform;
		if (TRANSFORM){setXLabel("Wave Vector (k)");}
		else {setXLabel("x");}
	}

	
}