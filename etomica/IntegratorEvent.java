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
public class IntegratorEvent {

    // Typed constants used to indicate the type of event integrator is
    // announcing
    public static final Type START =      new Type("Start",      1); //simulation is starting
    public static final Type INITIALIZE = new Type("Initialize", 2); //integrator is initializing
    public static final Type INTERVAL =   new Type("Interval",   4); //routine interval event
    public static final Type DONE =       new Type("Done",       8); //simulation is finished
    
    private final Type type;
    private Integrator source;

    public IntegratorEvent(Integrator source, Type t) {
        this.source = source;
        type = t;
    }

    public Type type() {
        return type;
    }
    
    public Integrator getSource() {
        return source;
    }

    //class used to mark the different types of interval events
    public final static class Type extends Constants.TypedConstant {
        public final int mask;
        private Type(String label, int mask) {
            super(label);
            this.mask = mask;
        }

        public static final Type[] choices = new Type[] {
                START, INTERVAL, DONE, INITIALIZE };

        public final Constants.TypedConstant[] choices() {
            return choices;
        }
    }
}
