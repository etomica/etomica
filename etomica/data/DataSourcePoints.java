package etomica.data;

import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataFunction;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataFunction.DataInfoFunction;
import etomica.data.types.DataPoint.DataInfoPoint;
import etomica.units.Dimension;
import etomica.units.Null;
import etomica.util.Function;

public class DataSourcePoints implements DataSource {

    private static final long serialVersionUID = 1L;
    private DataTag tag;
	private IDataInfo dataInfo;
    private DataDoubleArray dataX = null;
    private DataDoubleArray dataY = null;

    public DataSourcePoints(String label, Dimension dimension) {

        tag = new DataTag();
        dataInfo = new DataInfoPoint(label, dimension, new int[] {2, 2});
        dataInfo.addTag(tag);
    }

    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    public DataTag getTag() {
        return tag;
    }

    /**
     * This needs to be implemented.  The DataPump calls this function
     * to get the data and passes it downstream.  So, this needs to
     * return the entire data set(x and y)?
     */
    public Data getData() {
    	double [] d = new double[dataX.getLength() + dataY.getLength()];
    	double[] x = dataX.getData();
    	double[] y = dataY.getData();
    	for(int idx = 0; idx < dataX.getLength(); idx++) {
    		d[idx] = x[idx];
    	}
    	for(int idx = 0; idx < dataX.getLength(); idx++) {
    		d[idx + dataX.getLength()] = y[idx];
    	}

    	DataDoubleArray data = new DataDoubleArray(new int[] {dataX.getLength() + dataY.getLength()}, d);
        return data;
    }

    public void update(double[] xpts, double[] ypts) {
        updatePts(xpts, ypts);
    }

    private void updatePts(double[] xpts, double[] ypts) {
    	if(xpts == null || ypts == null || xpts.length != ypts.length) {
    		dataX = null;
    		dataY = null;
    		dataX = new DataDoubleArray(0);
    		dataY = new DataDoubleArray(0);
    	}
    	else {
	    	dataX = new DataDoubleArray(xpts.length/* * pointDim*/);
	    	dataY = new DataDoubleArray(ypts.length);
            double[] x = dataX.getData();
            double[] y = dataX.getData();

	    	for(int row = 0; row < xpts.length; row++) {
	    	    x[row] = xpts[row];
	    	    y[row] = ypts[row];
	    	}
    	}
    }

}
