/*
 * History
 * Created on Jul 28, 2004 by kofke
 */
package etomica;

import etomica.units.Dimension;

/**
 * @author kofke
 *
 * Accumulator that simply counts the number of times its 
 * <code>add</code> method is invoked.
 */
public class AccumulatorCounter extends Accumulator {

	/**
	 * @param parentElement
	 * @param dataSource
	 */
	public AccumulatorCounter(SimulationElement parentElement,
			DataSource dataSource) {
		super(parentElement, dataSource);
	}

	/* (non-Javadoc)
	 * @see etomica.Accumulator#add(double[])
	 */
	public void add(double[] values) {
		count++;

	}

	/* (non-Javadoc)
	 * @see etomica.DataSource#getData()
	 */
	public double[] getData() {
		data[0] = count;
		return data;
	}

	public void setNData(int nData) {
		// TODO Auto-generated method stub

	}

	/**
	 * Returns QUANTITY dimension.
	 */
	public Dimension getDimension() {
		return Dimension.QUANTITY;
	}

	int count;
	double[] data = new double[1];
}
