package etomica.graphics;
import etomica.DataSink;
import etomica.DataSource;
import etomica.EtomicaElement;
import etomica.EtomicaInfo;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageSegment;
import etomica.data.AccumulatorHistory;
import etomica.data.DataSourceUniform;
import etomica.data.meter.MeterPressureHard;
import etomica.simulations.HSMD2D;
import etomica.units.Dimension;
import etomica.units.Unit;
import etomica.utility.Arrays;

/**
 * Class for creating a plot of simulation data.
 *
 * @author David Kofke
 */

/* History
 * 01/01/03 (DAK) Modified doUpdate to plot only if isShowing
 * 05/18/04 (DAK) Added setSize method
 */
 
public class DisplayPlot extends Display implements EtomicaElement {
    
    private Plot plot;
    private javax.swing.JPanel panel = new javax.swing.JPanel();
    private boolean doLegend = true;
    protected DataGroup[] data = new DataGroup[0];
    protected DataGroup x;
    private DataSource xSource;
   
    public DisplayPlot() {
        super();
        plot = new Plot();
        panel.add(plot);
        xSource = new DataSourceUniform();
        x = new DataGroup(xSource.getData().length);
        setName("Data Plot");
        Thread runner = new Thread() {
            public void run() {
                while(true) {
                    doUpdate();
                    try {
                        Thread.sleep(100);
                    } catch(InterruptedException ex) {}
                }
            }
        };
        runner.start();
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("X-Y graphical plot of data");
        return info;
    }
    
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
     * Performs actions appropriate to addition or change of data source.
     * Implementation of abstract method from parent class.
     */
    public void setupDisplay() {
        panel.remove(plot);
        plot = new Plot();
        panel.add(plot);
        panel.revalidate();
        panel.repaint();
        
        for(int i=0; i<data.length; i++) {
            plot.addLegend(i,doLegend ? data[i].label : "");
        }
        
        setLabel(data[0].label);
        plot.setYLabel(data[0].label);

        plot.setXLabel(x.label + " ("+x.unit.symbol()+")");
        
        doUpdate();
    }
            
    public void doUpdate() {
        if(!plot.isShowing()) return;
        int nSource = data.length;
        plot.clear(false);
        x.y = xSource.getData();
        for(int k=0; k<nSource; k++) {
            for(int i=0; i<data[k].y.length; i++) {
                plot.addPoint(k, x.unit.fromSim(x.y[i]), data[k].unit.fromSim(data[k].y[i]), true);
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
     * Extend superclass method to update label with change of unit.
     */
    public void setXUnit(Unit u) {
        x.unit = u;
        if(plot != null && x.unit != null) plot.setXLabel(x.label + " ("+x.unit.symbol()+")");
    }
    
    public void setXSource(DataSource xSource) {
        this.xSource = xSource;
        x.y = xSource.getData();
    }
    
    public DataSink makeDataSink(Unit u) {
        DataGroup newGroup = (DataGroup)makeDataSink();
        newGroup.unit = u;
        return newGroup;
    }
    
    public DataSink makeDataSink() {
        DataGroup newGroup = new DataGroup(x.y.length);
        data = (DataGroup[])Arrays.addObject(data, newGroup);
        return newGroup;
    }
    
    public void removeDataSink(DataSink sink) {
        data = (DataGroup[])Arrays.removeObject(data, sink);
    }
    
    public void removeAllSinks() {
        data = new DataGroup[0];
    }

    public void setSize(int width, int height) {
        plot.setSize(width, height);
        panel.setSize(width, height);
    }

 /**
  * Define inner class as extension of ptolemy.plot.Plot
  * Does not override anything, but may want to later
  */
    public class Plot extends ptolemy.plot.Plot {
        
        public Plot() {
        	super();
        	setOpaque(false);
        }
        public void setTopPadding(int i) {
            _topPadding = i;
        }
    }
    
    private class DataGroup implements DataSink {
        double[] y;
        Dimension dimension = Dimension.UNDEFINED;
        Unit unit = Unit.NULL;
        String label = "";
        
        DataGroup(int n) {
            y = new double[n];
        }

        public synchronized void putData(double[] values) {
            if(y.length != values.length) y = (double[])values.clone();
            else System.arraycopy(values, 0, y, 0, values.length);
 //           doUpdate();
        }

        public void setDimension(Dimension dimension) {
            this.dimension = dimension;
        }

        public void setLabel(String label) {
            this.label = label;
        }
    }

    public static void main(String[] args) {
        HSMD2D sim = new HSMD2D();
        SimulationGraphic graphic = new SimulationGraphic(sim);
        sim.integrator.setIsothermal(true);
        AccumulatorHistory history = new AccumulatorHistory();
        history.setHistoryLength(1000);
        MeterPressureHard pressureMeter = new MeterPressureHard(sim.integrator);
        pressureMeter.setPhase(sim.phase);
        AccumulatorAverageSegment segment = new AccumulatorAverageSegment(
                pressureMeter, sim.integrator, 
                new AccumulatorAverage.Type[] {AccumulatorAverage.AVERAGE},
                new DisplayBox());
        segment.getDataPump().addDataSink(history);
        DisplayPlot plot = new DisplayPlot();
        history.addDataSink(plot.makeDataSink());
        plot.setXSource(history.getXSource());
        graphic.add(plot);

        graphic.makeAndDisplayFrame();
    }
}//end of class