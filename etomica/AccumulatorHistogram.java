/*
 * History
 * Created on Aug 4, 2004 by kofke
 */
package etomica;

import etomica.units.Dimension;
import etomica.utility.Histogram;

/**
 * @author kofke
 *
 */
public class AccumulatorHistogram extends Accumulator {
	
	Histogram[] histogram = new Histogram[0];
	int nData, nDataMinus1;
	/**
	 * 
	 */
	public AccumulatorHistogram() {
		super(Dimension.NULL);
		setNData(0);
	}


	/* (non-Javadoc)
	 * @see etomica.Accumulator#add(double[])
	 */
	public void add(double[] values) {
		if(values.length != nData) setNData(values.length);
		for(int i=nDataMinus1; i>=0; i--) {       		
			histogram[i].addValue(values[i]);
		}
	}

	/* (non-Javadoc)
	 * @see etomica.DataSource#getData()
	 */
	public double[] getData() {
		// TODO Auto-generated method stub
		return null;
	}

	/* (non-Javadoc)
	 * @see etomica.DataSource#getTranslator()
	 */
	public DataTranslator getTranslator() {
		// TODO Auto-generated method stub
		return null;
	}
	
	private void setNData(int nData) {
    	this.nData = nData;
    	nDataMinus1 = nData-1;
    	sum = redimension(nData, sum); 

}
