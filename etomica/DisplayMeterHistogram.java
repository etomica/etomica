//includes a main method to demonstrate and test use
package simulate;

import simulate.units.*;
import simulate.utility.Histogram;
import ptolemy.plot.Plot;
import javax.swing.JButton;
import java.awt.event.*;

/**
 * Class for creating a plot of histogram data kept by a meter.
 */
public class DisplayMeterHistogram extends Display implements Meter.User, MeterFunction.User {
    
    Plot plot;  //inner class defined below
    
    Meter meter;  //one or the other of these is used depending on meter type
    MeterFunction meterFunction;
    
    Unit xUnit, yUnit;
    private JButton resetButton = new JButton("Reset histogram");
    private boolean makeResetButton = false;
    private boolean multiHistogram = false;
   
    public DisplayMeterHistogram() {
        this(Simulation.instance);
    }
    public DisplayMeterHistogram(Simulation sim) {
        super(sim);
        setXUnit(new Unit(BaseUnit.Null.UNIT));
        setYUnit(new Unit(BaseUnit.Null.UNIT));
        setLabel("Plot");
    }
    
    public boolean isMakeResetButton() {return makeResetButton;}
    public void setMakeResetButton(boolean b) {
        if(b && !makeResetButton) { //asking for resetButton and not already made
            add(resetButton);
        }
        else if(!b) remove(resetButton);
        makeResetButton = b;
    }
    
    public Meter getMeter() {return meter;}
    public MeterFunction getMeterFunction() {return meterFunction;}
    public void setMeter(Meter m) {
        m.setHistogramming(true);
        if(plot != null) remove(plot);
        meter = m;
        meterFunction = null;
        multiHistogram = false;
        plot = new Plot();
        add(plot);
        setXUnit(m.defaultIOUnit());
        plot.setXLabel(meter.getLabel()+" ("+xUnit.symbol()+")");
        plot.setYLabel("Frequency");
        setLabel(meter.getLabel()); //tabbed-pane text
//        resetButton.addActionListener(new ActionListener() {
//            public void actionPerformed(ActionEvent evt) {meter.reset();}
//        });
    }
    public void setMeterFunction(MeterFunction m) {
        m.setHistogramming(true);
        if(plot != null) remove(plot);
        meterFunction = (MeterFunction)m;
        meter = null;
        multiHistogram = true;
        plot = new Plot();
        add(plot);
        setXUnit(m.defaultIOUnit());
        plot.setXLabel(meter.getLabel()+" ("+xUnit.symbol()+")");
        plot.setYLabel("Frequency");
        setLabel(meter.getLabel()); //tabbed-pane text
//        resetButton.addActionListener(new ActionListener() {
//            public void actionPerformed(ActionEvent evt) {meter.reset();}
//        });
    }
    
    public void setXUnit(Unit u) {xUnit = u;}
    public Unit getXUnit() {return xUnit;}
    public void setYUnit(Unit u) {yUnit = u;}
    public Unit getYUnit() {return yUnit;}
    
    public void doUpdate() {
        if(meter == null) return;
        plot.update();
    }
    
    public class Plot extends ptolemy.plot.Plot {
        
        public Plot() {
            setTitle("Histogram");
//            setYRange(0, 1);
//            setXRange(xMin, xMax);
            setImpulses(true); 
            setMarksStyle("none");
            update();
        }
        
        public void update() {  //redraws plot with current values of histogram
            clear(false);
            repaint();
//            setXRange(xMin, xMax);
            if(multiHistogram) {
                Histogram[] histograms = meterFunction.histogram();
                for(int i=0; i<histograms.length; i++) {
                    double[] hist = histograms[i].getHistogram();
                    double[] x = histograms[i].getX();
                    for(int j=0; j<hist.length; j++) {
                        addPoint(i,x[j],hist[j],false);
                    }
                }
            }
            else {
                Histogram histogram = meter.histogram();
                double[] hist = histogram.getHistogram();
                double[] x = histogram.getX();
                for(int j=0; j<hist.length; j++) {
                    addPoint(0,x[j],hist[j],false);
                }
            }
            repaint();
        }
    }//end of Plot inner class
    
    /**
     * Demonstrates how this class is implemented.
     */
    public static void main(String[] args) {
        java.awt.Frame f = new java.awt.Frame();   //create a window
        f.setSize(600,350);
        Simulation.makeSimpleSimulation();  //for more general simulations, replace this call with
                                            //construction of the desired pieces of the simulation
        //part that is unique to this demonstration
        Meter meter = new MeterPressureHard();
        Phase phase = Simulation.phase(0);
        DisplayMeterHistogram hist = new DisplayMeterHistogram();
        //end of unique part
                                            
		Simulation.instance.elementCoordinator.go(); //invoke this method only after all elements are in place
		                                    //calling it a second time has no effect
		                                    
        meter.setPhase(phase);//must come after go() because phase needs to have integrator for this call
        hist.setMeter(meter);
        hist.setUpdateInterval(10);
        f.add(Simulation.instance);         //access the static instance of the simulation to
                                            //display the graphical components
        f.pack();
        f.show();
        f.addWindowListener(new java.awt.event.WindowAdapter() {   //anonymous class to handle window closing
            public void windowClosing(java.awt.event.WindowEvent e) {System.exit(0);}
        });
    }

    
    
}