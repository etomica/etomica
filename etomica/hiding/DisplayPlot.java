package simulate;

import simulate.units.*;
import ptolemy.plot.Plot;

/**
 * Class for creating a plot of simulation data.
 * Data are obtained from a MeterFunction class that is identified via the
 * setMeter method.
 */
public class DisplayPlot extends Display {
    
    Plot plot;
    MeterFunction meter;
    Unit xUnit, yUnit;
    boolean useCurrentValue;  //flag indicating if display should show meter's current value (if true) or average (if false: default value)
   
    public DisplayPlot() {
        super();
        plot = new Plot();
        add(plot);
        setXUnit(Unit.Null.UNIT);
        setYUnit(Unit.Null.UNIT);
        setUseCurrentValue(false);
    }
    
    public void setMeter(MeterFunction m) {
        meter = m;
        plot.setXLabel(meter.getXLabel()+" ("+xUnit.symbol()+")");
        plot.setYLabel(meter.getLabel());
        setUseCurrentValue(getUseCurrentValue());
        if(m instanceof MeterMultiFunction) {
            int nF = ((MeterMultiFunction)m).nFunctions();
            for(int f=0; f<nF; f++) {
                MeterFunction m1 = ((MeterMultiFunction)m).meter(f);
                plot.addLegend(f,m1.getLabel());
            }
        }
    }
    
    public void setXUnit(Unit u) {xUnit = u;}
    public Unit getXUnit() {return xUnit;}
    public void setYUnit(Unit u) {yUnit = u;}
    public Unit getYUnit() {return yUnit;}
    
    /**
     * Specifies whether lines join points in plot (true) or not (false)
     */
    public void setConnected(boolean b) {plot.setConnected(b);}
    public boolean getConnected() {return plot.getConnected();}
    
    public void doUpdate() {
        plot.clear(false);
//        plot.repaint();
        if(meter instanceof MeterMultiFunction) {
            MeterMultiFunction m = (MeterMultiFunction)meter;
            for(int f=0; f<m.nFunctions(); f++) {
                double[] x = m.X(f);
                double[] y = (useCurrentValue) ? m.currentValue(f) : m.average(f);
                for(int i=0; i<x.length; i++) {
                    plot.addPoint(f,xUnit.fromSim(x[i]),yUnit.fromSim(y[i]),true);
                }
            }
        }
        else {
            double[] x = meter.X();
            double[] y = (useCurrentValue) ? meter.currentValue() : meter.average();
            for(int i=0; i<x.length; i++) {
                plot.addPoint(0,xUnit.fromSim(x[i]),yUnit.fromSim(y[i]),true);
            }
        }
        plot.repaint();
    }
    
    /**
     * Sets flag indicating if plot should be of instantaneous (current) value or running average
     */
    public void setUseCurrentValue(boolean b) {
        useCurrentValue = b;
        if(meter == null) return;
//        if(!useCurrentValue && !meter.isActive()) {meter.setActive(true);}
        if(!useCurrentValue && !meter.isActive()) {System.out.println("Warning: setting to use averages but meter is not active");}
    }
    public boolean getUseCurrentValue() {return useCurrentValue;}
 
 /**
  * Define inner class as extension of ptolemy.plot.Plot
  * Does not override anything, but may want to later
  */
    public class Plot extends ptolemy.plot.Plot {
    }
}