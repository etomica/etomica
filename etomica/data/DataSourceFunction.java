package etomica.data;

import java.io.Serializable;

import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataFunction;
import etomica.units.Dimension;
import etomica.utility.Function;

/**
 * Datasource formed as a wrapper around a function.  Uses a DataSourceUniform
 * to generate x values, and takes a Function that computes
 * the y values from them.  This DataSource returns a DataGroup with
 * two Data components; the first (0) is x, and the second (1) is y.
 * Useful for displaying a fixed function on a plot.
 */
 
 /* History
  * 09/08/02 (DAK) new
  */
  
public class DataSourceFunction implements DataSource, Serializable {
    
    public DataSourceFunction() {
        this(new Function.Constant(0.0));
    }
    public DataSourceFunction(Function function) {
        this("y(x)", Dimension.NULL, function, 100);
    }
    
    public DataSourceFunction(String label, Dimension dimension, Function function, int nValues) {
        xSource = new DataSourceUniform("x",Dimension.NULL,nValues,0,1);
        yData = new DataDoubleArray(label, dimension, nValues);
        this.function = function;
        update();
    }
    
    public DataInfo getDataInfo() {
        return data.getDataInfo();
    }
    
    /**
     * Returns the DataSourceUniform instance that generates the x values.
     * The range and spacing of the x values can be adjusted through the
     * methods of the returned instance.
     */
    public DataSourceUniform getXSource() {
        return xSource;
    }

    /**
     * Returns the DataFunction made by this source.
     */
    public Data getData() { 
        return data;
    }
    
    /**
     * @return Returns the function.
     */
    public Function getFunction() {
        return function;
    }
    /**
     * @param function The function to set.
     */
    public void setFunction(Function function) {
        this.function = function;
        update();
    }
    /**
     * Recalculates the y values from the current x values.  This must be
     * invoked if the methods of the DataSourceUniform given by getXSource
     * are used to change the x values.
     *
     */
    public void update() {
        boolean needUpdate = false;
        if (xData != (DataDoubleArray)xSource.getData()) {
            xData = (DataDoubleArray)xSource.getData();
            needUpdate = true;
        }
        double[] x = xData.getData();
        if (yData.getLength() != x.length) {
            yData = new DataDoubleArray(yData.getDataInfo(), new int[] {x.length});
            needUpdate = true;
        }
        if (needUpdate) {
            data = new DataFunction(new DataDoubleArray[] {xData}, yData);
        }
        double[] y = yData.getData();
        for(int i=0; i<x.length; i++) {
            y[i] = function.f(x[i]);
        }
    }
    
    private DataFunction data;
    private final DataSourceUniform xSource;
    private DataDoubleArray xData, yData;
    private Function function;
}//end of DataSourceFunction
