/*
 * History
 * Created on Jul 28, 2004 by kofke
 */
package etomica.data;

import etomica.Accumulator;
import etomica.DataTranslator;
import etomica.units.Dimension;

/**
 * Accumulator that simply counts the number of times its 
 * <code>add</code> method is invoked.
 */
public class AccumulatorCounter extends Accumulator {

	/**
	 * @param parentElement
	 * @param dataSource
	 */
	public AccumulatorCounter() {
		super(Dimension.QUANTITY);
	}

	/* (non-Javadoc)
	 * @see etomica.Accumulator#add(double[])
	 */
	public void putData(double[] values) {
		count++;
	}

	/* (non-Javadoc)
	 * @see etomica.DataSource#getData()
	 */
	public double[] getData() {
		data[0] = count;
		return data;
	}
	
	public void reset() {
		count = 0;
	}

	int count;
	final double[] data = new double[1];
	
	public DataTranslator getTranslator() {return DataTranslator.IDENTITY;}
}
