package etomica.graphics;
import etomica.*;
import java.awt.*;
import java.awt.event.*;
import java.util.Observer;
import java.util.Observable;
import java.beans.Beans;
import java.beans.PropertyChangeSupport;
import java.beans.PropertyChangeListener;
import javax.swing.JPanel;

/**
 * Superclass of all classes that display something from the simulation.  
 * Included are displays of graphical and tabular data, and views of the molecules as they
 * move about during the simulation.  Displays are also used to output simulation data and
 * results to file.
 *
 * @author David Kofke
 */
public abstract class Display extends SimulationElement implements GraphicalElement, Integrator.IntervalListener, java.io.Serializable {

    public static final String VERSION = "Display:01.07.25";
    
  /**
   * Number of interval events received between updates of display.
   * Default value is 1.
   */
    int updateInterval;
  /**
   * Counter used to track number of interval events since last update of display.
   */
    protected int iieCount;
  /**
   * Descriptive text for the display.  Used to set the text in the tab of the tabbed display panel.
   */
   protected String label;
   
   private JPanel panel = new JPanel();
    
    protected Phase phase;  //consider removing this and putting in subclasses only as needed
                            //used at least by DisplayPhase, DisplayTable, DisplayScrollingGraph, DisplayToConsole
    
    // Constructor
    public Display(Simulation sim) {
        super(sim, Display.class);
	    setUpdateInterval(1);
    }
    
    /**
     * Method to associate a phase with this display.
     * Default action is to simply set the phase field to the given value
     */
    public void setPhase(Phase p) {phase = p;}  //2D needed to manipulate dimensions array directly
    /**
     * @return the phase associated with this display
     */
    public final Phase phase() {return phase;}
    
    /**
     * Accessor method for phase.  Equivalent to phase() method.
     */
    public final Phase getPhase() {return phase();}
    /**
     * Method of Simulation.GraphicElement interface.
     * Default action is to return this Display (which is a Panel) as the graphic object.
     * May override in subclass to return a more appropriate graphical element, or none at all.
     */
    public Component graphic(Object obj) {return panel;}
    
    /**
     * Same as graphic method with a null argument.
     */
    public final Component graphic() {return graphic(null);}
    
    public void initialize() {}
    
    /**
     * Method called to update the display.  
     * This method is called after the display receives 
     * an integrator IntervalEvent updateInterval times.
     * After completing this method the Display does a repaint.
     */
    public abstract void doUpdate();
    
    public void repaint() {panel.repaint();}
    
    /**
     * Method of Integrator.IntervalListener interface.
     * After receiving an event updateInterval times, the method will invoke the 
     * doUpdate method and then call repaint().
     */
    public void intervalAction(Integrator.IntervalEvent evt) {
	    if(evt.type() == Integrator.IntervalEvent.INITIALIZE || --iieCount == 0) {
	        iieCount = updateInterval;
	        doUpdate();
            repaint();
	    }
    }

    /**
     * Accessor method for the update interval
     * Number of interval events received between updates of display.
     */
    public final int getUpdateInterval() {return updateInterval;}
    /**
     * Accessor method for the update interval
     * Number of interval events received between updates of display.
     */
    public final void setUpdateInterval(int i) {
        if(i > 0) {
            updateInterval = i;
            iieCount = updateInterval+1;
        }
    }
    
    /**
     * Accessor method of the label describing the display.
     * 
     * @return The given label
     */
    public String getLabel() {return label;}

    /**
     * Accessor method of the label describing the display.
     * 
     * @param label The label string describing the display
     */
    public void setLabel(String label) {
        String oldLabel = this.label;
        this.label = label;
        support.firePropertyChange("label", oldLabel, label);
    }
    /**
     * Method to set the name of this object
     * 
     * @param name The name string to be associated with this object
     */
    public final void setName(String name) {
        String oldName = this.name;
        super.setName(name);
        this.name = name;
        support.firePropertyChange("name", oldName, name);
    }

    public void addPropertyChangeListener(PropertyChangeListener listener) {
        support.addPropertyChangeListener(listener);
    }
    public void removePropertyChangeListener(PropertyChangeListener listener) {
        support.removePropertyChangeListener(listener);
    }
    
    protected PropertyChangeSupport support = new PropertyChangeSupport(this);
}
