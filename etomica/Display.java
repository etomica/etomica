package etomica;
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
 * results to file
 */
public abstract class Display extends Panel implements Simulation.GraphicalElement, Integrator.IntervalListener, java.io.Serializable {

    public static final String VERSION = "Display:01.03.11.0";
    private String name;
    
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
    
    protected Phase phase;  //consider removing this and putting in subclasses only as needed
                            //used at least by DisplayPhase, DisplayTable, DisplayScrollingGraph, DisplayToConsole
    
    private final Simulation parentSimulation;
    private boolean added = false;
    
    // Constructor
    public Display(Simulation sim) {
        parentSimulation = sim;
	    setUpdateInterval(1);
        parentSimulation.register(this);
    }
    
    public final Simulation parentSimulation() {return parentSimulation;}
    public final Class baseClass() {return Display.class;}
    public final boolean wasAdded() {return added;}
    public final void setAdded(boolean b) {added = b;}
    
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
    public Component graphic(Object obj) {return this;}
    
    /**
     * Method called to update the display.  
     * This method is called after the display receives 
     * an integrator IntervalEvent updateInterval times.
     * After completing this method the Display does a repaint.
     */
    public abstract void doUpdate();
    
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
     * Accessor method of the name of this object
     * 
     * @return The given name
     */
    public final String getName() {return name;}

    /**
     * Method to set the name of this object
     * 
     * @param name The name string to be associated with this object
     */
    public final void setName(String name) {
        String oldName = this.name;
        this.name = name;
        support.firePropertyChange("name", oldName, name);
    }

    /**
     * Overrides the Object class toString method to have it return the output of getName
     * 
     * @return The name given to the object
     */
    public String toString() {return getName();}  //override Object method
          
    
    public void addPropertyChangeListener(PropertyChangeListener listener) {
        support.addPropertyChangeListener(listener);
    }
    public void removePropertyChangeListener(PropertyChangeListener listener) {
        support.removePropertyChangeListener(listener);
    }
    
    protected PropertyChangeSupport support = new PropertyChangeSupport(this);
}
