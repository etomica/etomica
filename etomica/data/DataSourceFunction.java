package etomica.data;

import etomica.Data;
import etomica.DataInfo;
import etomica.DataSource;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataFunction;
import etomica.units.Dimension;
import etomica.utility.Function;
import etomica.utility.IntegerRange;

/**
 * Datasource formed as a wrapper around a function.  Takes a
 * data source for the x values, and a Function that computes
 * the y values from them.  Default x source is DataSourceUniform.
 * Useful for displaying a fixed function on a plot.
 */
 
 /* History
  * 09/08/02 (DAK) new
  */
  
public class DataSourceFunction implements DataSource, DataSourceDependent, java.io.Serializable {
    
    public DataSourceFunction() {
        this(new Function.Constant(0.0));
    }
    public DataSourceFunction(Function function) {
        data = new DataFunction(new DataInfo("y",Dimension.NULL),new DataInfo("x",Dimension.NULL));
        xSource = new DataSourceUniform();
        this.function = function;
        update();
    }
    
    public DataInfo getDataInfo() {
        return data.getDataInfo();
    }
    
    public DataSource getXSource() {return xSource;}
    
    public DataSource[] getDataSource() {
    	return new DataSource[] {xSource};
    }
    public void setDataSource(DataSource[] source) {
    	xSource = source[0];
    }
    public IntegerRange dataSourceCountRange() {
    	return new IntegerRange(1,1);
    }
    
    public Data getData() {
        return data;
    }
    
    public void update() {
        DataDoubleArray xData = (DataDoubleArray)xSource.getData();
        data.setLength(xData.getLength());
        data.getTData().E(xData);
        y = data.getData();
        for(int i=0; i<x.length; i++) {
            y[i] = function.f(x[i]);
        }
    }
    
    private DataFunction data;
    private DataSource xSource;
    private Function function;
    private double[] x, y;
}//end of DataSourceFunction
