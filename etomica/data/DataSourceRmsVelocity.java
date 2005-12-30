package etomica.data;

import etomica.EtomicaElement;
import etomica.atom.Atom;
import etomica.atom.iterator.AtomIterator;
import etomica.space.ICoordinateKinetic;
import etomica.units.Undefined;

/**
 * Meter for the root-mean-square velocity of a set of atoms. Useful to obtain
 * histograms of atom speeds.
 * 
 * @author David Kofke
 */

public class DataSourceRmsVelocity extends DataSourceScalar implements EtomicaElement {

	//TODO define a dimension for this property
	public DataSourceRmsVelocity() {
		super("RMS Velocity", Undefined.DIMENSION);
	}

	/**
	 * Returns the rms velocity of the atoms given by the iterator. Value is
	 * given in the first element of the array, which is always of dimension 1.
	 */
	public double getDataAsScalar() {
		iterator.reset();
		int count = 0;
		double value = 0.0;
		while (iterator.hasNext()) {
			Atom atom = iterator.nextAtom();
			value += Math.sqrt(((ICoordinateKinetic)atom.coord).velocity().squared());
			count++;
		}
		value /= count;
		return value;
	}
    
    public int getDataLength() {
        return 1;
    }

	/**
	 * Sets the iterator defining the atoms for which the RMS velocity is
	 * calculated.
	 * 
	 * @param iterator
	 */
	public void setIterator(AtomIterator iterator) {
		if (iterator == null)
			this.iterator = AtomIterator.NULL;
		else
			this.iterator = iterator;
	}

	/**
	 * @return The iterator defining the atoms used for the RMS velocity
	 *         calculation
	 */
	public AtomIterator getIterator() {
		return iterator;
	}

	private AtomIterator iterator = AtomIterator.NULL;

}//end of DataSourceVelocityRms
