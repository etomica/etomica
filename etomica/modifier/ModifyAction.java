package etomica.modifier;

import etomica.Action;
import etomica.Modifier;
import etomica.units.Dimension;


/**
 * Wraps a Modifier instance so that it can be passed to the Controller
 * (or other activity group) for execution there.  Permits a modification
 * initiated by a user interface thread to be executed by the controller
 * thread (for example).
 *
 * @author David Kofke
 *
 */

/*
 * History
 * Created on Feb 1, 2005 by kofke
 */
public class ModifyAction implements Action, Modifier {

    /**
     * Constuctor requires the wrapped Modifier, which is final.
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

    /* (non-Javadoc)
     * @see etomica.Modifier#setValue(double)
     */
    public void setValue(double newValue) {
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

    private final Modifier modifier;
    private double value;
}
