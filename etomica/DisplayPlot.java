package etomica;

import etomica.units.*;
import ptolemy.plot.Plot;
import javax.swing.JButton;
import java.awt.event.*;

/**
 * Class for creating a plot of simulation data.
 * Data are obtained from a DataSource that is identified via the
 * setDataSource method.
 *
 * @author David Kofke
 */
 
 //need to update to move away from meter and more toward data source
 
public class DisplayPlot extends Display implements /*DataSource.User,*/ EtomicaElement {
    
    public String getVersion() {return "DisplayPlot:01.05.23/"+Display.VERSION;}

    private Plot plot;
    private DataSource[] ySource;
    private DataSource xSource;
    private int nSource;
    private Unit xUnit, yUnit;
    private DataSource.ValueType whichValueX, whichValue;
    private javax.swing.JPanel panel = new javax.swing.JPanel();
   
    public DisplayPlot() {
        this(Simulation.instance);
    }
    public DisplayPlot(Simulation sim) {
        super(sim);
        plot = new Plot();
        panel.add(plot);
        setXUnit(Unit.NULL);
        setYUnit(Unit.NULL);
        setName("Data Plot");
//        new Thread(this).start();
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("X-Y graphical plot of data");
        return info;
    }
    
    /**
     * Accessor method to plot class so that its properties can be edited.
     */
    public final Plot getPlot() {return plot;}
        
    /**
     * Overrides superclass method to return a JPanel holding the plot.
     */
    public java.awt.Component graphic(Object obj) {return panel;}

    /**
     * @deprecated.  Use setDataSource instead.
     */
    public void setMeterFunction(MeterFunction m) {setDataSource(m);}
    
    public DataSource[] getDataSource() {return ySource;}
    
    public DataSource getDataSource(int i) {
        return (ySource != null && i < nSource) ? ySource[i] : null;
    }
    
    public void setDataSource(DataSource s) {
        setDataSource(new DataSource[] {s});
    }
    
    public void setDataSource(DataSource[] s) {
        ySource = s;
        if(s == null) {nSource = 0; return;}
        nSource = s.length;
        if(nSource == 0) return;
        for(int i=0; i<nSource; i++) {
            plot.addLegend(i,ySource[i].getLabel());
        }
        setLabel(ySource[0].getLabel());
        plot.setYLabel(ySource[0].getLabel());
        //change unit if dimension of new source is different from current source        
        if(yUnit.dimension() != ySource[0].getDimension()) 
            setYUnit(ySource[0].getDimension().defaultIOUnit());
    }
    
    public void setXSource(DataSource s) {
        xSource = s;
        if(xSource == null) {
            setXUnit(Unit.NULL);
            plot.setXLabel("");
        }
        else {
            if(xUnit.dimension() != xSource.getDimension())
                setXUnit(xSource.getDimension().defaultIOUnit());
            plot.setXLabel(xSource.getLabel()+" ("+xUnit.symbol()+")");
        }
    }
        
    public void setXUnit(Unit u) {xUnit = u;}
    public Unit getXUnit() {return xUnit;}
    public void setYUnit(Unit u) {yUnit = u;}
    public Unit getYUnit() {return yUnit;}
    
    public void doUpdate() {
        if(ySource == null) return;
        plot.clear(false);
        if(xSource != null) {
            for(int k=0; k<nSource; k++) {
                double[] x = xSource.values(null);
                double[] y = ySource[k].values(whichValue);
                for(int i=0; i<y.length; i++) {
                    plot.addPoint(k,xUnit.fromSim(x[i]),yUnit.fromSim(y[i]),true);
                }//for i
            }//for k
        }//end if
        else {//xSource is null
            for(int k=0; k<nSource; k++) {
                double[] y = ySource[k].values(whichValue);
                for(int i=0; i<y.length; i++) {
                    plot.addPoint(k,(double)i,yUnit.fromSim(y[i]),true);
                }//for i
            }//for k
        }//end else
        plot.repaint();
    }//end doUpdate method
    
    /**
     * Sets whether meter displays average, current value, last block average, etc.
     */
    public void setWhichValue(DataSource.ValueType type) {
        whichValue = type;
    }
    public DataSource.ValueType getWhichValue() {return whichValue;}
     
 /**
  * Define inner class as extension of ptolemy.plot.Plot
  * Does not override anything, but may want to later
  */
    public class Plot extends ptolemy.plot.Plot {
    }
}