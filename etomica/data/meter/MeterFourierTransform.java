/*
 * Created on Jul 7, 2003
 */
package etomica.data.meter;


import etomica.Atom;
import etomica.Default;
import etomica.Phase;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.atom.iterator.AtomIteratorPhaseDependent;
import etomica.units.Dimension;
import etomica.utility.*;


/**
 * @author Mitchell Tai
 *
 * Meter which monitors the component k values of a function as a result of
 * a Fourier Transform. 
 * 
 * What you say!
 */
public class MeterFourierTransform extends MeterFunction {
	AtomIteratorPhaseDependent iterator; 
	int nMolecules;
	double[] data;
	public boolean REAL=true;		// toggle to meter the real/imaginary transform
	public boolean TRANSFORM=true;	// toggle to transform
	FastFourierTransform fourier;
	

	public MeterFourierTransform(boolean real) {
		super();
		this.REAL=real;
		if (real){setLabel("Real");}
		else {setLabel("Imaginary");}
		xDataSource.setLabel("Wave Vector (k)");
		fourier = new FastFourierTransform();
	}
	public void setPhase(Phase[] p) {
		super.setPhase(p);
		iterator = new AtomIteratorLeafAtoms();
		setNDataPerPhase(p[0].atomCount());	
		nMolecules = nDataPerPhase;
	}

	public Dimension getXDimension() {
		return Dimension.NULL;
	}
	/**
	 * returns the array of real or imaginary numbers
	 * 
	 * No doubt 
	 */
	public double[] getDataAsArray(Phase aPhase) {
        iterator.setPhase(aPhase);
		iterator.reset();
		Atom currentAtom;
		data = new double[nMolecules];
		for (int i=0;i<nMolecules;i++) {
			currentAtom = iterator.nextAtom();
			// calculate the displacement from their home position.
			data[i] = currentAtom.coord.position().x(0) - (i*Default.BOX_SIZE/nMolecules);
			
			// accomodates for excessive distances due to wrap around
			if (data[i]>Default.BOX_SIZE/2){data[i]-=Default.BOX_SIZE;}
			if (data[i]<-Default.BOX_SIZE/2){data[i]+=Default.BOX_SIZE;}
		}
		fourier.setData(data);
		if(TRANSFORM)fourier.transform();
//		x=fourier.getIndex();
		if (REAL){return fourier.getReal();}
		return fourier.getImaginary();
	}

	public Dimension getDimension() {
		return Dimension.NULL;
	}

	public boolean isTRANSFORM() {
		return TRANSFORM;
	}

	public void setTRANSFORM(boolean transform) {
		TRANSFORM = transform;
		if (TRANSFORM){xDataSource.setLabel("Wave Vector (k)");}
		else {xDataSource.setLabel("x");}
	}

	
}