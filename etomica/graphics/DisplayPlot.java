package etomica.graphics;
import etomica.*;
import etomica.units.*;

/**
 * Class for creating a plot of simulation data.
 *
 * @author David Kofke
 */

/* History
 * 01/01/03 (DAK) Modified doUpdate to plot only if isShowing
 * 05/18/04 (DAK) Added setSize method
 */
 
public class DisplayPlot extends DisplayDataSources implements EtomicaElement {
    
    private Plot plot;
    private javax.swing.JPanel panel = new javax.swing.JPanel();
    private boolean doLegend = true;
   
    public DisplayPlot() {
        this(Simulation.instance);
    }
    public DisplayPlot(Simulation sim) {
        super(sim);
        plot = new Plot();
        panel.add(plot);
        setName("Data Plot");
        setWhichValue(MeterAbstract.AVERAGE);
//        new Thread(this).start();
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
        if(ySource == null) return;
        if(ySource.length > 1) {
            for(int i=0; i<ySource.length; i++) plot.addLegend(i,doLegend ? ySource[i].getLabel() : "");
        }
        setLabel(ySource[0].getLabel());
        plot.setYLabel(ySource[0].getLabel());

        if(xSource == null) {
            plot.setXLabel("");
        }
        else {
            if(xUnit.dimension() != xSource.getDimension())
            plot.setXLabel(getXLabel() + " ("+xUnit.symbol()+")");
        }
        doUpdate();
    }
            
    public void doUpdate() {
    	if(!plot.isShowing()) return;
        super.doUpdate();
        if(ySource == null) return;
        int nSource = ySource.length;
        plot.clear(false);
        ///new stuff
        for(int k=0; k<nSource; k++) {
            for(int i=0; i<x.length; i++) {
              if(!Double.isNaN(y[k][i])) { 
				plot.addPoint(k, xUnit.fromSim(x[i]), yUnit[k].fromSim(y[k][i]), true);
              } else if(i==x.length-1) {
				plot.addPoint(k, xUnit.fromSim(x[i]), 0.0, false);
              }
            }//for i
        }//for k
        plot.repaint();
        ///end of new stuff
/*        if(xSource != null) {
            for(int k=0; k<nSource; k++) {
                double[] x = xSource.values(whichValueX);
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
       
        plot.repaint(); */
    }//end doUpdate method
    
    /**
     * Extend superclass method to update label with change of unit.
     */
    public void setXUnit(Unit u) {
        xUnit = u;
        if(plot != null && xUnit != null) plot.setXLabel(getXLabel() + " ("+xUnit.symbol()+")");
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

}//end of class