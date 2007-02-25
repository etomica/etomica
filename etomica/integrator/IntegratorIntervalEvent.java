package etomica.integrator;

public class IntegratorIntervalEvent extends IntegratorEvent {

    private static final long serialVersionUID = 1L;
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
        private static final long serialVersionUID = 1L;

        protected IntervalEventType(String label) {
            super(label);
        }
        
        public static Type[] choices() {
            return new Type[] {INTERVAL};
        }
    }

    
}