package etomica.integrator;

import etomica.util.Constants;
import etomica.util.EnumeratedType;

/**
 * Event object given to registered listeners when an integrator fires
 * any type event.  Object provides information about the type of event
 * being fired, and gives a reference to the integrator.
 *
 * @author David Kofke
 *
 */

/*
 * History
 * Created on Apr 27, 2005 by kofke
 */
public class IntegratorEvent implements java.io.Serializable {

    // Typed constants used to indicate the type of event integrator is
    // announcing
    public static final Type START =      new Type("Start",      1); //simulation is starting
    public static final Type INITIALIZE = new Type("Initialize", 2); //integrator is initializing
    public static final Type INTERVAL =   new Type("Interval",   4); //routine interval event
    public static final Type DONE =       new Type("Done",       8); //simulation is finished
    
    private final Type type;
    private Integrator source;

    /**
     * 
     */
    IntegratorEvent(Integrator source, Type type) {
        this.source = source;
        this.type = type;
    }

    public Type type() {
        return type;
    }
    
    public Integrator getSource() {
        return source;
    }

    //class used to mark the different types of interval events
    public final static class Type extends EnumeratedType {
        public final int mask;
        private Type(String label, int mask) {
            super(label);
            this.mask = mask;
        }

        public static final Type[] choices = new Type[] {
                START, INTERVAL, DONE, INITIALIZE };

        public final EnumeratedType[] choices() {
            return choices;
        }
    }

}
