/*
 * History
 * Created on Jul 28, 2004 by kofke
 */
package etomica.data;

import etomica.DataTranslator;
import etomica.units.Dimension;

/**
 * Data processor that simply counts the number of times its 
 * <code>putData</code> method is invoked.
 */
public class AccumulatorCounter extends DataAccumulator {

	/**
	 * @param parentElement
	 * @param dataSource
	 */
	public AccumulatorCounter() {
		setDimension(Dimension.QUANTITY);
	}

	/* (non-Javadoc)
	 * @see etomica.Accumulator#add(double[])
	 */
	protected void addData(double[] values) {
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
