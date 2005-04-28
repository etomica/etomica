package etomica;

/**
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 *
 * @author David Kofke
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
