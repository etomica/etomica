package etomica.graphics;
import ptolemy.plot.Plot;
import etomica.EtomicaElement;
import etomica.EtomicaInfo;
import etomica.data.DataSinkTable;
import etomica.data.DataSource;
import etomica.data.DataSourceUniform;
import etomica.data.DataTableListener;
import etomica.data.types.DataArithmetic;
import etomica.data.types.DataTable;
import etomica.units.Null;
import etomica.units.Unit;
import etomica.units.systems.UnitSystem;
import etomica.util.Arrays;

/**
 * Class for creating a plot of simulation data.  Data for plot is taken
 * from a DataSinkTable instance, which can be obtained via the getDataTable
 * method.  Data sets are added to the plot by piping the data to a new DataSink
 * obtained from the makeColumn method of the data table.  
 *
 * @author David Kofke
 */

public class DisplayPlot extends Display implements DataTableListener, EtomicaElement {
    
    /**
     * Creates a plot with a new, empty, DataSinkTable.
     */
    public DisplayPlot() {
        this(new DataSinkTable());
    }
    
    /**
     * Creates a plot using data from the given table.
     */
    public DisplayPlot(DataSinkTable table) {
        super();
        this.dataTable = table;
        table.addTableListener(this);
        plot = new Plot();
        panel.add(plot);
        x = new DataSourceUniform();
        setName("Data Plot");
        units = new Unit[table.getColumnCount()];
        for(int i=0; i<units.length; i++) {
            units[i] = table.getColumn(i).getDimension().getUnit(UnitSystem.SIM);
        }
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("X-Y graphical plot of data");
        return info;
    }
    
    /**
     * Returns the DataSinkTable instance that specifies the plot data.
     * Data sets are added to the plot by piping them to a new data sink
     * given by the makeColumn method of the table.
     * @return
     */
    public DataSinkTable getDataTable() {
        return dataTable;
    }
    
    /**
     * Causes the display of the plot to be updated.
     */
    public void tableDataChanged(DataSinkTable table) {
        doUpdate();
    }
    
    /**
     * Updates the units array for the new column, using the default units.
     */
    public void tableColumnCountChanged(DataSinkTable table) {
        int oldColumnCount = units.length;
        int newColumnCount = table.getColumnCount();
        if(newColumnCount > oldColumnCount) {
            units = (Unit[])Arrays.resizeArray(units, newColumnCount);
            for(int i=oldColumnCount; i<newColumnCount; i++) {
                units[i] = (defaultUnit != null) ? defaultUnit : table.getColumn(i).getDimension().getUnit(UnitSystem.SIM);
            }
        } else {
            //TODO have DisplayTable adjust appropriately to removal of columns; this works only for newColumnCount = 0
            //may need to remember columns and reassign units to match
            units = (Unit[])Arrays.resizeArray(units, newColumnCount);
        }
    }

    /**
     * Has no effect. Part of the DataTableListener interface.
     */
    public void tableRowCountChanged(DataSinkTable table) {
        //do nothing
    }
    
    /**
     * Redraws the plot.
     */
    public void doUpdate() {
        if(!plot.isShowing()) return;
        int nSource = dataTable.getColumnCount();
        plot.clear(false);
        DataArithmetic xValues = (DataArithmetic)x.getData();
        for(int k=0; k<nSource; k++) {
            final double[] data = dataTable.getColumn(k).getData();
            for(int i=0; i<data.length; i++) {
                double y = units[k].fromSim(data[i]);
                if (!Double.isNaN(y)) {
                    plot.addPoint(k, xUnit.fromSim(xValues.getValue(i)), units[k].fromSim(data[i]), true);
                }
            }
        }
        plot.repaint();

    }//end doUpdate method

    /**
     * Mutator for flag determining if a legend is to be shown.
     * Default is true.
     */
    public void setDoLegend(boolean b) {
        doLegend = b;
        for(int i=0; i<dataTable.getColumnCount(); i++) {
            plot.removeLegend(i);
            plot.addLegend(i,b ? dataTable.getColumn(i).getHeading(): "");
        }
    }
    
    /**
     * Accessor for flag determining if a legend is to be shown.
     * Default is true.
     */
    public boolean isDoLegend() {return doLegend;}
    
    /**
     * Accessor method to plot class so that its properties can be edited.
     */
    public final Plot getPlot() {return plot;}
        
    /**
     * Overrides superclass method to return a JPanel holding the plot.
     */
    public java.awt.Component graphic(Object obj) {return panel;}


    /**
     * Extend superclass method to update label with change of unit.
     */
    public void setXUnit(Unit u) {
        xUnit = u;
        if(plot != null && xUnit != null) plot.setXLabel(x.getDataInfo().getLabel() + " ("+xUnit.symbol()+")");
    }
    
    /**
     * Sets the source for the x-axis value.  Default is a DataSourceUniform instance.
     */
    public void setXSource(DataSource xSource) {
        x = xSource;
    }
    
    /**
     * Returns the source for the x-axis value.
     */
    public DataSource getXSource() {
        return x;
    }
    
    /**
     * Sets the display units of the given column to the given unit.
     */
    public void setUnit(DataTable.Column column, Unit newUnit) {
        for(int i=0; i<dataTable.getColumnCount(); i++) {
            if(dataTable.getColumn(i) == column) {
                units[i] = newUnit;
                break;
            }
        }
    }
    
    public void setUnit(Unit newUnit) {
    	defaultUnit = newUnit;
        for(int i=0; i<dataTable.getColumnCount(); i++) {
            units[i] = newUnit;
        }
    }

    /**
     * Sets the drawn size of the plot
     * @param width new width, in pixels
     * @param height new height, in pixels
     */
    public void setSize(int width, int height) {
        plot.setSize(width, height);
        panel.setSize(width, height);
    }
    
    private final DataSinkTable dataTable;
    private Plot plot;
    private javax.swing.JPanel panel = new javax.swing.JPanel();
    private boolean doLegend = true;
    protected DataSource x;
    private Unit xUnit = Null.UNIT;
    private Unit[] units = new Unit[0];
    private Unit defaultUnit;
}