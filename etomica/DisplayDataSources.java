package etomica;

import java.awt.Component;
import javax.swing.JTable;
import javax.swing.table.*;
import javax.swing.Box;
import javax.swing.JScrollPane;
import etomica.units.Unit;

/**
 * Parent class for displays that show information from one or more DataSource objects.
 * Examples are DisplayPlot and DisplayTableFunction.  This class implements the methods
 * of the DataSource.MultiUser interface, which enable setting and changing the 
 * data source objects.
 *
 * @author David Kofke
 */
public abstract class DisplayDataSources extends Display implements DataSource.MultiUser
{
    public String getVersion() {return "DisplayDataSources:01.05.29/"+Display.VERSION;}

    protected DataSource[] ySource;
    protected DataSource xSource;
    protected Unit yUnit, xUnit;
    protected double x[];
    protected double[][] y;
    protected DataSource.ValueType whichValueX, whichValue;
        
    public DisplayDataSources(Simulation sim)  {
        super(sim);
        setXUnit(Unit.NULL);
        setYUnit(Unit.NULL);
    }
    
    /**
     * Method called whenever a data source is added or removed.
     */
    protected abstract void setupDisplay();
    
    /**
     * Returns the array of all data sources being tabulated.
     */
    public DataSource[] getDataSources() {return ySource;}
    
    /**
     * Sets the ith data source in the array of data sources.  Takes
     * no action if i is greater than the current number of data sources (minus 1).
     */
    public DataSource getDataSources(int i) {
        if(ySource == null) return null;
        if(i < 0 || i >= ySource.length) throw new ArrayIndexOutOfBoundsException();
        else return ySource[i];
    }
    /**
     * Sets the array of data sources being tabulated.  Existing sources
     * are discarded.
     */
    public void setDataSources(DataSource[] s) {
        if(s == null || s.length == 0 || s[0] == null) {
            ySource = null;
            return;
        }
        ySource = s;
        y = new double[ySource.length][];
        //change unit if dimension of new source is different from current source        
        if(yUnit.dimension() != ySource[0].getDimension()) 
            setYUnit(ySource[0].getDimension().defaultIOUnit());
        if(xSource == null && ySource[0] instanceof DataSource.X) 
            xSource = ySource[0];
        setupX();
        setupDisplay();
    }
    
    public void setDataSources(int i, DataSource s) {
        if(ySource == null && i != 0) throw new NullPointerException();
        else if(ySource == null && i == 0) {setDataSources(s); return;}
        if(i < 0 || i >= ySource.length) throw new ArrayIndexOutOfBoundsException();
        else ySource[i] = s;
        setupX();
        setupDisplay();
    }
    
    /**
     * Sets the given source as the only data source being tabulated.
     * Existing sources are discarded.
     */
    public void setDataSources(DataSource s) {
        setDataSources(new DataSource[] {s});
    }
    
    /**
     * Adds the given source to the sources being plotted.
     * Existing sources are retained.
     */
    public void addDataSources(DataSource s) {
        if(s == null) return;
        int nSource = (ySource == null) ? 0 : ySource.length;
        DataSource[] newSources = new DataSource[nSource+1];
        for(int i=0; i<nSource; i++) newSources[i] = ySource[i];
        newSources[nSource] = s;
        setDataSources(newSources);
    }
    
    /**
     * Adds the given sources to the sources being plotted.
     * Existing sources are retained.
     */
    public void addDataSources(DataSource[] s) {
        if(s == null) return;
        int nSource = (ySource == null) ? 0 : ySource.length;
        DataSource[] newSources = new DataSource[nSource+s.length];
        for(int i=0; i<nSource; i++) newSources[i] = ySource[i];
        for(int i=nSource; i<newSources.length; i++) newSources[i] = s[i-nSource];
        setDataSources(newSources);
    }
    
    public void setXSource(DataSource s) {
        xSource = s;
        if(xSource == null) setXUnit(Unit.NULL);
        setupX();
        setupDisplay();
    }
    
    /**
     * Sets up array of x values.
     */
    private void setupX() {
        if(xSource != null) {
            if(ySource[0] instanceof DataSource.X) {
                x = ((DataSource.X)ySource[0]).xValues();
 //               setXLabel(((DataSource.X)xSource).getXLabel());
                setXUnit(((DataSource.X)ySource[0]).getXDimension().defaultIOUnit());
            }
            else
                x = xSource.values(whichValueX);                
        }
        else { //no xSource specified; take sequentially
            x = new double[maxYPoints()];
            for(int i=0; i<x.length; i++) {x[i] = i;}
        }
    }
    
    /**
     * Returns the number of points for the data source having the most points.
     */
    protected int maxYPoints() {
        if(ySource == null || ySource.length == 0) return 0;
        else {
            int k = 0;
            for(int i=0; i<ySource.length; i++) {
                double[] ys = ySource[i].values(whichValue);
                int n = (ys != null) ? ys.length : 0;
                k = (k < n) ? n : k;
            }
            return k;
        }
    }
    
    /**
     * Sets variable indicating which value is to be taken from data sources.
     */
    public void setWhichValue(DataSource.ValueType type) {
        whichValue = type;
        setupDisplay();
    }
    public DataSource.ValueType getWhichValue() {return whichValue;}
     
    /**
     * Sets which value to use from x data source.
     */
    public void setWhichValueX(DataSource.ValueType type) {
        whichValueX = type;
        setupDisplay();
    }
    public DataSource.ValueType getWhichValueX() {return whichValueX;}
    
    public void setXUnit(Unit u) {xUnit = u;}
    public Unit getXUnit() {return xUnit;}
    public void setYUnit(Unit u) {yUnit = u;}
    public Unit getYUnit() {return yUnit;}
    
    public void doUpdate() {
        //update all y values at once so not to invoke calculation of all
        //y values when just one of them is accessed
        for(int i=0; i<ySource.length; i++) {
           y[i] = ySource[i].values(whichValue);
        }
    }
}//end of DisplayDataSources