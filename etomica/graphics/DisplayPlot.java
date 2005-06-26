package etomica.graphics;
import etomica.DataSource;
import etomica.EtomicaElement;
import etomica.EtomicaInfo;
import etomica.data.DataBin;
import etomica.data.DataSourceUniform;
import etomica.data.DataTable;
import etomica.data.DataTableListener;
import etomica.data.types.DataArithmetic;
import etomica.units.Unit;
import etomica.utility.Arrays;

/**
 * Class for creating a plot of simulation data.  Data for plot is taken
 * from a DataTable instance, which can be obtained via the getDataTable
 * method.  Data sets are added to the plot by piping the data to a new DataSink
 * obtained from the makeColumn method of the data table.  
 *
 * @author David Kofke
 */

public class DisplayPlot extends Display implements DataTableListener, EtomicaElement {
    
    /**
     * Creates a plot with a new, empty, DataTable.
     */
    public DisplayPlot() {
        this(new DataTable());
    }
    
    /**
     * Creates a plot using data from the given table.
     */
    public DisplayPlot(DataTable table) {
        super();
        this.dataTable = table;
        table.addTableListener(this);
        plot = new Plot();
        panel.add(plot);
        x = new DataSourceUniform();
        setName("Data Plot");
        units = new Unit[table.getColumnCount()];
        for(int i=0; i<units.length; i++) {
            units[i] = table.getColumn(i).getDataInfo().getDimension().defaultIOUnit();
        }
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("X-Y graphical plot of data");
        return info;
    }
    
    /**
     * Returns the DataTable instance that specifies the plot data.
     * Data sets are added to the plot by piping them to a new data sink
     * given by the makeColumn method of the table.
     * @return
     */
    public DataTable getDataTable() {
        return dataTable;
    }
    
    /**
     * Causes the display of the plot to be updated.
     */
    public void tableDataChanged(DataTable table) {
        doUpdate();
    }
    
    /**
     * Updates the units array for the new column, using the default units.
     */
    public void tableColumnAdded(DataTable table, DataBin newColumn) {
        units = (Unit[])Arrays.addObject(units, newColumn.getDataInfo().getDimension().defaultIOUnit());
    }

    /**
     * Causes the corresponding units element to be removed.
     */
    public void tableColumnRemoved(DataTable table, int index, DataBin oldColumn) {
        units = (Unit[])Arrays.removeObject(units, units[index]);
    }
    
    /**
     * Has no effect. Part of the DataTableListener interface.
     */
    public void tableRowCountChanged(DataTable table, int oldCount, int newCount) {
        //do nothing
    }
    
    /**
     * Performs actions appropriate to addition or change of data source.
     * Implementation of abstract method from parent class.
     */
    public void setupDisplay() {
        panel.remove(plot);
        plot = new Plot();
        panel.add(plot);
        panel.revalidate();
        panel.repaint();
        
        for(int i=0; i<dataTable.getColumnCount(); i++) {
            plot.addLegend(i,doLegend ? dataTable.getColumn(i).getDataInfo().getLabel(): "");
        }
        
        setLabel(x.getDataInfo().getLabel());
        if(dataTable.getColumnCount() > 0) {
            plot.setYLabel(dataTable.getColumn(0).getDataInfo().getLabel());
        }

        plot.setXLabel(x.getDataInfo().getLabel()+ " ("+xUnit.symbol()+")");
        
        doUpdate();
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
            DataArithmetic column = (DataArithmetic)dataTable.getColumn(k).getData();
            for(int i=0; i<column.getLength(); i++) {
                plot.addPoint(k, xUnit.fromSim(xValues.getValue(i)), units[k].fromSim(column.getValue(i)), true);
//              if(!Double.isNaN(data[k].y[i])) { 
//                plot.addPoint(k, x.unit.fromSim(x.y[i]), data[k].unit.fromSim(data[k].y[i]), true);
//              } else if(i==x.y.length-1) {
//				plot.addPoint(k, x.unit.fromSim(x.y[i]), 0.0, false);
//              }
            }//for i
        }//for k
        plot.repaint();

    }//end doUpdate method

    /**
     * Mutator for flag determining if a legend is to be shown.
     * Default is true.
     */
    public void setDoLegend(boolean b) {
        doLegend = b;
        setupDisplay();
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
     * Sets the display units of the given column to the given unit.
     */
    public void setUnit(DataBin column, Unit newUnit) {
        for(int i=0; i<dataTable.getColumnCount(); i++) {
            if(dataTable.getColumn(i) == column) {
                units[i] = newUnit;
                break;
            }
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
    
    private final DataTable dataTable;
    private Plot plot;
    private javax.swing.JPanel panel = new javax.swing.JPanel();
    private boolean doLegend = true;
    protected DataSource x;
    private Unit xUnit = Unit.NULL;
    private Unit[] units = new Unit[0];

 /**
  * Define inner class as extension of ptolemy.plot.Plot
  * Does not override anything, but may want to later
  */
    public class Plot extends ptolemy.plot.Plot {
        
        public Plot() {
        	super();
//        	setOpaque(false);
        }
        public void setTopPadding(int i) {
            _topPadding = i;
        }
    }
    
//    public static void main(String[] args) {
//        HSMD2D sim = new HSMD2D();
//        SimulationGraphic graphic = new SimulationGraphic(sim);
//        sim.integrator.setIsothermal(true);
//        AccumulatorHistory history = new AccumulatorHistory();
//        history.setHistoryLength(1000);
//        MeterPressureHard pressureMeter = new MeterPressureHard(sim.space,sim.integrator);
//        pressureMeter.setPhase(sim.phase);
//        AccumulatorAverageSegment segment = new AccumulatorAverageSegment(
//                pressureMeter, sim.integrator, 
//                new AccumulatorAverage.Type[] {AccumulatorAverage.AVERAGE},
//                new DisplayBox());
//        segment.getDataPump().addDataSink(history);
//        DisplayPlot plot = new DisplayPlot();
//        history.addDataSink(plot.getDataTable().makeColumn());
//        plot.setXSource(history.getXSource());
//        graphic.add(plot);
//
//        graphic.makeAndDisplayFrame();
//    }
}//end of class