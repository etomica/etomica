package etomica.data;

import java.io.Serializable;

import etomica.Data;
import etomica.DataInfo;
import etomica.DataSource;
import etomica.data.types.DataDoubleArray;
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
        yData = new DataDoubleArray("y", Dimension.NULL);
        xSource = new DataSourceUniform();
        xData = (DataDoubleArray)xSource.getData();
        xData.getDataInfo().setLabel("x");
        data = new DataGroup("y(x)", Dimension.NULL, new Data[] {xData,yData});
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
     * Returns a 2-element DataGroup with the x values in getData(0), and
     * the y values in getData(1).
     */
    public Data getData() { 
        return data;
    }
    
    /**
     * Recalculates the y values from the current x values.  This must be
     * invoked if the methods of the DataSourceUniform given by getXSource
     * are used to change the x values.
     *
     */
    public void update() {
        double[] x = xData.getData();
        yData.setLength(x.length);
        double[] y = yData.getData();
        for(int i=0; i<x.length; i++) {
            y[i] = function.f(x[i]);
        }
    }
    
    private final DataGroup data;
    private final DataSourceUniform xSource;
    private final DataDoubleArray xData, yData;
    private Function function;
}//end of DataSourceFunction
