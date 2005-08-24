package etomica.modifier;

import etomica.action.Action;
import etomica.units.Dimension;


/**
 * Wraps a Modifier instance so that the modification can be performed
 * as the consequence of an action. Permits for example a modification
 * initiated by a user-interface thread to be executed by the controller
 * thread.  Call to setValueForAction stores value for modifier, which
 * is then invoked when a subsequent call to actionPerformed is made. 
 *
 * @author David Kofke
 *
 */

/*
 * History
 * Created on Feb 1, 2005 by kofke
 */
public class ModifyAction implements Action, java.io.Serializable {

    /**
     * Constuctor requires the wrapped Modifier, which is final.
     * Default activity group is null, and default targetAction is this.
     */
    public ModifyAction(Modifier modifier) {
        this.modifier = modifier;
    }

    /**
     * Executes the setValue method of the wrapped modifier, using
     * the value given in the most recent call to this.setValue. 
     * Uses zero if this.setValue was not previously invoked.
     */
    public void actionPerformed() {
        modifier.setValue(value);
    }

    /* (non-Javadoc)
     * @see etomica.Modifier#getLabel()
     */
    public String getLabel() {
        return modifier.getLabel();
    }

    /**
     * Stores the new value to be set by the wrapped modifier.  Value is 
     * not set until actionPerformed is invoked.
     */
    public void setValueForAction(double newValue) {
        value = newValue;
    }

    /* (non-Javadoc)
     * @see etomica.Modifier#getValue()
     */
    public double getValue() {
        return modifier.getValue();
    }

    /* (non-Javadoc)
     * @see etomica.Modifier#getDimension()
     */
    public Dimension getDimension() {
        return modifier.getDimension();
    }
    
    public Modifier getWrappedModifier() {
        return modifier;
    }
        
    private final Modifier modifier;
    private double value;
}
