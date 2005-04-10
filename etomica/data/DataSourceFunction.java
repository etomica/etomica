package etomica.data;

import etomica.DataSource;
import etomica.DataTranslator;
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
  
public class DataSourceFunction implements DataSource, DataSourceDependent {
    
    public DataSourceFunction() {
        this(new Function.Constant(0.0));
    }
    public DataSourceFunction(Function function) {
        xSource = new DataSourceUniform();
        this.function = function;
        dimension = Dimension.NULL;
        label = "";
        update();
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
    
    public double[] getData() {
        return y;
    }
    
    public int getDataLength() {
        return y.length;
    }
    
    public double[] xValues() {return x;}
    
    public void update() {
        x = xSource.getData();
        y = new double[x.length];
        for(int i=0; i<x.length; i++) {
            y[i] = function.f(x[i]);
        }
    }
    
    public void setDimension(Dimension dimension) {this.dimension = dimension;}
    public Dimension getDimension() {return dimension;}
    
    public Dimension getXDimension() {return xSource.getDimension();}
    
    public String getXLabel() {return xSource.getLabel();}
    public String getLabel() {return label;}
    public void setLabel(String label) {this.label = label;}
    
    private DataSource xSource;
    private Function function;
    private double[] x, y;
    private String label;
    private Dimension dimension;
    
    public DataTranslator getTranslator() {return DataTranslator.IDENTITY;}
 
}//end of DataSourceFunction