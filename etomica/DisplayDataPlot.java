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
public class DisplayDataPlot extends Display implements DataSource.User, EtomicaElement {
    
    public String getVersion() {return "DisplayDataPlot:01.03.24/"+Display.VERSION;}

    Plot plot;
    DataSource ySource, xSource;
    Unit xUnit, yUnit;
    DataSource.ValueType whichValueX, whichValue;
    private JButton resetButton = new JButton("Reset averages");
    private boolean makeResetButton = false;
   
    public DisplayDataPlot() {
        this(Simulation.instance);
    }
    public DisplayDataPlot(Simulation sim) {
        super(sim);
        plot = new Plot();
        add(plot);
        setXUnit(new Unit(BaseUnit.Null.UNIT));
        setYUnit(new Unit(BaseUnit.Null.UNIT));
        setName("Data Plot");
//        setWhichValueX(MeterAbstract.ValueType.AVERAGE);
        setWhichValue(MeterAbstract.ValueType.AVERAGE);
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("X-Y graphical plot of data");
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
    
    /**
     * @deprecated.  Use setDataSource instead.
     */
    public void setMeterFunction(MeterFunction m) {setDataSource(m);}
    
    public DataSource getDataSource() {return ySource;}
    
    public void setDataSource(DataSource s) {
        ySource = s;
        if(ySource == null) return;
        if(ySource instanceof MeterFunction) {
            if(!(whichValue instanceof MeterAbstract.ValueType)) whichValue = MeterAbstract.ValueType.AVERAGE;
            final MeterFunction meter = (MeterFunction)ySource;
            plot.setXLabel(meter.getXLabel()+" ("+xUnit.symbol()+")");
            plot.setYLabel(meter.getLabel());
            setWhichValue(getWhichValue());
            if(ySource instanceof MeterMultiFunction) {
                int nF = ((MeterMultiFunction)meter).nFunctions();
                for(int k=0; k<nF; k++) {
                    MeterFunction m1 = ((MeterMultiFunction)meter).meter(k);
                    plot.addLegend(k,m1.getLabel());
                }
            }
            setLabel(meter.getLabel()); //tabbed-pane text
            resetButton.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent evt) {meter.reset();}
            });
        }
        else if(ySource instanceof Meter) {
            if(!(whichValue instanceof Meter.ValueType)) whichValue = Meter.ValueType.HISTORY;
            final Meter meter = (Meter)ySource;
            plot.setYLabel(meter.getLabel());
            setWhichValue(getWhichValue());
            setLabel(meter.getLabel());
            resetButton.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent evt) {meter.reset();}
            });
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
        if(ySource == null) return;
        plot.clear(false);
//        plot.repaint();
        if(ySource instanceof MeterFunction) {
            MeterFunction meter = (MeterFunction)ySource;
            if(meter instanceof MeterMultiFunction) {
                MeterMultiFunction m = (MeterMultiFunction)meter;
                for(int k=0; k<m.nFunctions(); k++) {
                    double[] x = m.X(k);
                    double[] y = m.values(k,(MeterAbstract.ValueType)whichValue);
                    for(int i=0; i<x.length; i++) {
                        plot.addPoint(k,xUnit.fromSim(x[i]),yUnit.fromSim(y[i]),true);
                    }
                }
            }
            else {
                double[] x = meter.X();
        //     double[] y = (useCurrentValue) ? meter.currentValue() : meter.average();
                double[] y = meter.values(whichValue);
                for(int i=0; i<x.length; i++) {
                    plot.addPoint(0,xUnit.fromSim(x[i]),yUnit.fromSim(y[i]),true);
                }
            }
        }
        else if(ySource instanceof Meter) {
            Meter meter = (Meter)ySource;
            double[] y = meter.values(whichValue);
            for(int i=0; i<y.length; i++) {
                plot.addPoint(0,(double)i, yUnit.fromSim(y[i]),true);
            }
        }
        plot.repaint();
    }
    
    /**
     * Sets whether meter displays average, current value, last block average, etc.
     */
    public void setWhichValue(DataSource.ValueType type) {
        whichValue = type;
        if(ySource == null) return;
        if(ySource instanceof MeterFunction) {
            MeterFunction meter = (MeterFunction)ySource;
            if(!(whichValue==MeterAbstract.ValueType.MOST_RECENT) && !meter.isActive())
                {System.out.println("Warning: setting to use averages but meter is not active");}
        }
    }
    public DataSource.ValueType getWhichValue() {return whichValue;}
     
 /**
  * Define inner class as extension of ptolemy.plot.Plot
  * Does not override anything, but may want to later
  */
    public class Plot extends ptolemy.plot.Plot {
    }
}