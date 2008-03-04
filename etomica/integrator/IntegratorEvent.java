package etomica.integrator;

import etomica.api.IIntegrator;
import etomica.util.EnumeratedType;

/**
 * Event object given to registered listeners when an integrator fires
 * any type event.  Object provides information about the type of event
 * being fired, and gives a reference to the integrator.
 *
 * @author David Kofke
 *
 */
public class IntegratorEvent implements java.io.Serializable {

    // Typed constants used to indicate the type of event integrator is
    // announcing
    
    private static final long serialVersionUID = 1L;
    private final Type type;
    private IIntegrator source;

    /**
     * 
     */
    IntegratorEvent(IIntegrator source, Type type) {
        this.source = source;
        this.type = type;
    }

    public Type type() {
        return type;
    }
    
    public IIntegrator getSource() {
        return source;
    }

    //class used to mark the different types of interval events
    public static class Type extends EnumeratedType {
        private static final long serialVersionUID = 1L;

        protected Type(String label) {
            super(label);
        }

        public static Type[] choices() {
            return new Type[]{IntegratorNonintervalEvent.RESET};
        }
    }

}
