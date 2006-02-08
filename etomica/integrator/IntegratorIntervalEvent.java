package etomica.integrator;

import etomica.util.EnumeratedType;


/*
 * History
 * Created on Apr 26, 2005 by kofke
 */
public class IntegratorIntervalEvent extends IntegratorEvent {

    private final int interval;
    
    public IntegratorIntervalEvent(Integrator source, int interval) {
        super(source, INTERVAL);
        this.interval = interval;
    }

    public int getInterval() {
        return interval;
    }

    public static final IntervalEventType INTERVAL =  new IntervalEventType("Interval"); //routine interval event

    public static class IntervalEventType extends Type {
        protected IntervalEventType(String label) {
            super(label);
        }
        
        public static Type[] choices() {
            return new Type[] {INTERVAL};
        }
    }

    
}