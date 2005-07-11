/*
 * History
 * Created on Jul 28, 2004 by kofke
 */
package etomica.data;

import etomica.Data;
import etomica.DataInfo;
import etomica.data.types.DataInteger;
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
        data = new DataInteger("Counter", Dimension.QUANTITY);
	}
    
    public DataInfo getDataInfo() {
        return data.getDataInfo();
    }

	/* (non-Javadoc)
	 * @see etomica.Accumulator#add(double[])
	 */
	protected void addData(Data dummyData) {
		data.x++;
	}

	/* (non-Javadoc)
	 * @see etomica.DataSource#getData()
	 */
	public Data getData() {
		return data;
	}
	
	public void reset() {
		data.x = 0;
	}
	
	final DataInteger data;
}
