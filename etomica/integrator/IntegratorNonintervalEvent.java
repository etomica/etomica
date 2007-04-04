package etomica.integrator;


/**
 * Event thrown by integrator when it announces reaching special points in the 
 * simulation process, such as its beginning and end. 
 *
 * @author David Kofke and Andrew Schultz
 *
 */
public class IntegratorNonintervalEvent extends IntegratorEvent {

    private static final long serialVersionUID = 1L;

    public IntegratorNonintervalEvent(Integrator source, NonintervalEventType type) {
        super(source, type);
    }

    public static final NonintervalEventType INITIALIZE = new NonintervalEventType("Initialize"); //integrator is initializing
    public static final NonintervalEventType DONE =       new NonintervalEventType("Done");       //simulation is finished

    public static class NonintervalEventType extends Type {
        private static final long serialVersionUID = 1L;

        protected NonintervalEventType(String label) {
            super(label);
        }
        
        public static Type[] choices() {
            return new NonintervalEventType[] {INITIALIZE,DONE};
        }
    }
}
