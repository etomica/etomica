package etomica.simulation;

import etomica.phase.Phase;
import etomica.util.EnumeratedType;

public class SimulationEvent extends java.util.EventObject {
    
    public SimulationEvent(Object source) {
    	super(source);
    }
    
    public void setPhase(Phase p) {
        phase = p;
    }
    public Phase getPhase() {
        return phase;
    }
    
    public void setType(Type t) {type = t;}
    public Type type() {return type;}
    
    protected Type type;
    protected Phase phase;

    public static class Type extends EnumeratedType {

        /**
         * Required to guarantee singleton when deserializing.
         * @return the singleton INSTANCE
         */
        private Object readResolve() {
            for (int i=0; i<CHOICES.length; i++) {
                if (this.toString().equals(CHOICES[i].toString())) {
                    return CHOICES[i];
                }
            }
            throw new RuntimeException("unknown PhaseEvent type: "+this);
        }

        private Type(String label) {super(label);}
        public static final Type[] CHOICES = new Type[] {
            new Type("Phase added"),
            new Type("Phase removed")};
        public final EnumeratedType[] choices() {return CHOICES;}
    }
    public static final Type PHASE_ADDED =        Type.CHOICES[0];
    public static final Type PHASE_REMOVED =      Type.CHOICES[1];
}