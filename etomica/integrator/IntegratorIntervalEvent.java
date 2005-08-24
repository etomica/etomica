package etomica.integrator;

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


}