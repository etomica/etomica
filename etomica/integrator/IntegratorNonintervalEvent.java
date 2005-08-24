package etomica.integrator;

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

    public IntegratorNonintervalEvent(Integrator source, Type type) {
        super(source, type);
        if(type == INTERVAL) throw new IllegalArgumentException("Interval event should not be constructed as an IntegratorNonintervalEvent");
    }

}
