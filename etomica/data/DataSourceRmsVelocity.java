package etomica.data;

import etomica.Atom;
import etomica.AtomIterator;
import etomica.space.ICoordinateKinetic;

/**
 * Meter for the root-mean-square velocity of a set of atoms. Useful to obtain
 * histograms of atom speeds.
 * 
 * @author David Kofke
 */

public class DataSourceRmsVelocity extends DataSourceAdapter {

	//TODO define a dimension for this property
	public DataSourceRmsVelocity() {
		super(etomica.units.Dimension.UNDEFINED);
	}

	/**
	 * Returns the rms velocity of the atoms given by the iterator. Value is
	 * given in the first element of the array, which is always of dimension 1.
	 */
	public double[] getData() {
		iterator.reset();
		int count = 0;
		value[0] = 0.0;
		while (iterator.hasNext()) {
			Atom atom = iterator.nextAtom();
			value[0] = atom.type.rm()
					* Math.sqrt(((ICoordinateKinetic)atom.coord).momentum().squared());
			count++;
		}
		value[0] /= (double) count;
		return value;
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

	private final double[] value = new double[1];

}//end of DataSourceVelocityRms
