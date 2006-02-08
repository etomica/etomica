package etomica.integrator;

import etomica.integrator.IntegratorEvent.Type;

/**
 * Event thrown by integrator when it announces reaching special points in the 
 * simulation process, such as its beginning and end. 
 *
 * @author David Kofke and Andrew Schultz
 *
 */

/*
 * History
 * Created on Apr 26, 2005 by kofke
 */
public class IntegratorNonintervalEvent extends IntegratorEvent {

    public IntegratorNonintervalEvent(Integrator source, NonintervalEventType type) {
        super(source, type);
    }

    public static final NonintervalEventType START =      new NonintervalEventType("Start");      //simulation is starting
    public static final NonintervalEventType INITIALIZE = new NonintervalEventType("Initialize"); //integrator is initializing
    public static final NonintervalEventType DONE =       new NonintervalEventType("Done");       //simulation is finished

    public static class NonintervalEventType extends Type {
        protected NonintervalEventType(String label) {
            super(label);
        }
        
        public static Type[] choices() {
            return new NonintervalEventType[] {START,INITIALIZE,DONE};
        }
    }
}
