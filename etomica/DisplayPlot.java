package etomica;

import etomica.units.*;
import ptolemy.plot.Plot;
import javax.swing.JButton;
import java.awt.event.*;

/**
 * Class for creating a plot of simulation data.
 * Data are obtained from a MeterFunction class that is identified via the
 * setMeter method.
 *
 * @author David Kofke
 */
public class DisplayPlot extends Display implements MeterFunction.User, EtomicaElement {
    
    public String getVersion() {return "DisplayPlot:01.03.24.0/"+Display.VERSION;}

    Plot plot;
    MeterFunction meter;
    Unit xUnit, yUnit;
    MeterAbstract.ValueType whichValue;
    private JButton resetButton = new JButton("Reset averages");
    private boolean makeResetButton = false;
   
    public DisplayPlot() {
        this(Simulation.instance);
    }
    public DisplayPlot(Simulation sim) {
        super(sim);
        plot = new Plot();
        add(plot);
        setXUnit(new Unit(BaseUnit.Null.UNIT));
        setYUnit(new Unit(BaseUnit.Null.UNIT));
        setName("Plot");
        setWhichValue(MeterAbstract.ValueType.AVERAGE);
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Plot of data given by a meter function");
        return info;
    }
    
    public boolean isMakeResetButton() {return makeResetButton;}
    public void setMakeResetButton(boolean b) {
        if(b && !makeResetButton) { //asking for resetButton and not already made
            add(resetButton);
        }
        else if(!b) remove(resetButton);
        makeResetButton = b;
    }
    
    public void setMeterFunction(MeterFunction m) {
        meter = m;
        if(meter == null) return;
        plot.setXLabel(meter.getXLabel()+" ("+xUnit.symbol()+")");
        plot.setYLabel(meter.getLabel());
        setWhichValue(getWhichValue());
        if(m instanceof MeterMultiFunction) {
            int nF = ((MeterMultiFunction)m).nFunctions();
            for(int f=0; f<nF; f++) {
                MeterFunction m1 = ((MeterMultiFunction)m).meter(f);
                plot.addLegend(f,m1.getLabel());
            }
        }
        setLabel(meter.getLabel()); //tabbed-pane text
        resetButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {meter.reset();}
        });
    }
    public MeterFunction getMeterFunction() {return meter;}
    
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
        if(meter == null) return;
        plot.clear(false);
//        plot.repaint();
        if(meter instanceof MeterMultiFunction) {
            MeterMultiFunction m = (MeterMultiFunction)meter;
            for(int f=0; f<m.nFunctions(); f++) {
                double[] x = m.X(f);
                double[] y = m.value(f,whichValue);
                for(int i=0; i<x.length; i++) {
                    plot.addPoint(f,xUnit.fromSim(x[i]),yUnit.fromSim(y[i]),true);
                }
            }
        }
        else {
            double[] x = meter.X();
       //     double[] y = (useCurrentValue) ? meter.currentValue() : meter.average();
            double[] y = meter.value(whichValue);
            for(int i=0; i<x.length; i++) {
                plot.addPoint(0,xUnit.fromSim(x[i]),yUnit.fromSim(y[i]),true);
            }
        }
        plot.repaint();
    }
    
    /**
     * Sets whether meter displays average, current value, last block average, etc.
     */
    public void setWhichValue(MeterAbstract.ValueType type) {
        whichValue = type;
        if(meter == null) return;
        if(!(whichValue==MeterAbstract.ValueType.MOST_RECENT) && !meter.isActive())
            {System.out.println("Warning: setting to use averages but meter is not active");}
    }
    public MeterAbstract.ValueType getWhichValue() {return whichValue;}
     
 /**
  * Define inner class as extension of ptolemy.plot.Plot
  * Does not override anything, but may want to later
  */
    public class Plot extends ptolemy.plot.Plot {
    }
}