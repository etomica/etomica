package etomica.graphics;
import java.awt.Component;
import java.beans.PropertyChangeListener;
import java.beans.PropertyChangeSupport;

import javax.swing.JPanel;

import etomica.Action;
import etomica.Phase;
import etomica.utility.NameMaker;

/**
 * Superclass of all classes that display something from the simulation.  
 * Included are displays of graphical and tabular data, and views of the molecules as they
 * move about during the simulation.  Displays are also used to output simulation data and
 * results to file.
 *
 * @author David Kofke
 */
public abstract class Display implements Action, GraphicalElement, java.io.Serializable {

  /**
   * Descriptive text for the display.  Used to set the text in the tab of the tabbed display panel.
   */
   protected String label;
   
   private JPanel panel = new JPanel();
    
    protected Phase phase;  //consider removing this and putting in subclasses only as needed
                            //used at least by DisplayPhase, DisplayTable, DisplayScrollingGraph, DisplayToConsole
    // Constructor
    public Display() {
        setName(NameMaker.makeName(this.getClass()));
    }
    
    public void actionPerformed() {
        repaint();
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
    
    public void repaint() {panel.repaint();}
    
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
        this.name = name;
        support.firePropertyChange("name", oldName, name);
    }
    /**
     * Accessor method of the name of this phase
     * 
     * @return The given name of this phase
     */
    public final String getName() {return name;}

    /**
     * Overrides the Object class toString method to have it return the output of getName
     * 
     * @return The name given to the phase
     */
    public String toString() {return getName();}


    public void addPropertyChangeListener(PropertyChangeListener listener) {
        support.addPropertyChangeListener(listener);
    }
    public void removePropertyChangeListener(PropertyChangeListener listener) {
        support.removePropertyChangeListener(listener);
    }

    protected PropertyChangeSupport support = new PropertyChangeSupport(this);
    
    private String name;
}
