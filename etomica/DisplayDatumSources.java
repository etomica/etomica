package etomica;

import java.awt.Component;
import javax.swing.JTable;
import javax.swing.table.*;
import javax.swing.Box;
import javax.swing.JScrollPane;
import etomica.units.Unit;

/**
 * Parent class for displays that show information from one or more DatumSource objects.
 * This class implements the methods
 * of the DatumSource.MultiUser interface, which enable setting and changing the 
 * Datum source objects.
 *
 * @author David Kofke
 */
public abstract class DisplayDatumSources extends Display implements DatumSource.MultiUser
{
    public String getVersion() {return "DisplayDatumSources:01.05.29/"+Display.VERSION;}

    protected DatumSource[] ySource;
    protected Unit yUnit;//ignored for now
    protected double[][] y;
    protected String[] labels;
    protected DataSource.ValueType[] whichValues;//values shown for every source
        
    public DisplayDatumSources(Simulation sim)  {
        super(sim);
        setYUnit(Unit.NULL);
        whichValues = new DataSource.ValueType[1]; //default whichValues[0] is null
    }
    
    /**
     * Method called whenever a Datum source is added or removed.
     */
    protected abstract void setupDisplay();
    
    /**
     * Returns the array of all Datum sources being tabulated.
     */
    public DatumSource[] getDatumSources() {return ySource;}
    
    /**
     * Sets the ith Datum source in the array of Datum sources.  Takes
     * no action if i is greater than the current number of Datum sources (minus 1).
     */
    public DatumSource getDatumSources(int i) {
        if(ySource == null) return null;
        if(i < 0 || i >= ySource.length) throw new ArrayIndexOutOfBoundsException();
        else return ySource[i];
    }
    /**
     * Sets the array of Datum sources being tabulated.  Existing sources
     * are discarded.
     */
    public void setDatumSources(DatumSource[] s) {
        ySource = s;
        if(s == null || s.length == 0) return;
        y = new double[ySource.length][whichValues.length];
        //change unit if dimension of new source is different from current source        
        if(yUnit.dimension() != ySource[0].getDimension()) 
            setYUnit(ySource[0].getDimension().defaultIOUnit());
        setupLabels();
        setupDisplay();
    }
    
    public void setDatumSources(int i, DatumSource s) {
        if(ySource == null && i != 0) throw new NullPointerException();
        else if(ySource == null && i == 0) {setDatumSources(s); return;}
        if(i < 0 || i >= ySource.length) throw new ArrayIndexOutOfBoundsException();
        else ySource[i] = s;
        setupLabels();
        setupDisplay();
    }
    
    /**
     * Sets the given source as the only Datum source being tabulated.
     * Existing sources are discarded.
     */
    public void setDatumSources(DatumSource s) {
        setDatumSources(new DatumSource[] {s});
    }
    
    /**
     * Adds the given source to the sources being plotted.
     * Existing sources are retained.
     */
    public void addDatumSources(DatumSource s) {
        if(s == null) return;
        int nSource = (ySource == null) ? 0 : ySource.length;
        DatumSource[] newSources = new DatumSource[nSource+1];
        for(int i=0; i<nSource; i++) newSources[i] = ySource[i];
        newSources[nSource] = s;
        setDatumSources(newSources);
    }
    
    /**
     * Adds the given sources to the sources being plotted.
     * Existing sources are retained.
     */
    public void addDatumSources(DatumSource[] s) {
        if(s == null) return;
        int nSource = (ySource == null) ? 0 : ySource.length;
        DatumSource[] newSources = new DatumSource[nSource+s.length];
        for(int i=0; i<nSource; i++) newSources[i] = ySource[i];
        for(int i=nSource; i<newSources.length; i++) newSources[i] = s[i-nSource];
        setDatumSources(newSources);
    }
    
    /**
     * Mutator for the array of strings that are used to label the datum values.
     * If the given array of strings is not of the same length as the current
     * array of ySource, the method returns without performing any action.
     */
    public void setLabels(String[] text) {
        if(text.length != ySource.length) return;
        labels = text;
        setupLabels();
        setupDisplay();
    }
    public String[] getLabels() {
        return labels;
    }
    
    private void setupLabels() {
        if(labels == null || (labels.length != ySource.length) && ySource != null) {
            labels = new String[ySource.length];
            for(int i=0; i<ySource.length; i++) {
                if(ySource[i] != null) labels[i] = ySource[i].getLabel();
            }
        }
    }
    
    /**
     * Sets variable indicating which value is to be taken from Datum sources.
     */
    public void setWhichValues(DataSource.ValueType[] types) {
        whichValues = types;
        y = new double[ySource.length][whichValues.length];
        setupDisplay();
    }
    public void setWhichValues(DataSource.ValueType type) {
        setWhichValues(new DataSource.ValueType[] {type});
    }
    public DataSource.ValueType[] getWhichValues() {return whichValues;}
     
    public void setYUnit(Unit u) {yUnit = u;}
    public Unit getYUnits() {return yUnit;}
    
    public void doUpdate() {
        if(ySource == null) return;
        //update all y values at once so not to invoke calculation of all
        //y values when just one of them is accessed
        for(int i=0; i<ySource.length; i++) {
            if(y[i] == null) continue;
            for(int k=0; k<whichValues.length; k++) {
               y[i][k] = ySource[i].value(whichValues[k]);
            }
        }
    }
}//end of DisplayDatumSources