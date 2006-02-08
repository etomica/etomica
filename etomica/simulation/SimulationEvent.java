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

    public static final Type PHASE_ADDED =   new Type("Phase added");
    public static final Type PHASE_REMOVED = new Type("Phase removed");

    public static class Type extends EnumeratedType {

        protected Type(String label) {super(label);}

        public static Type[] choices() {
            return new Type[] {PHASE_ADDED,PHASE_REMOVED};
        }

        /**
         * Required to guarantee singleton when deserializing.
         * @return the singleton INSTANCE
         */
        private Object readResolve() {
            Type[] choices = choices();
            for (int i=0; i<choices.length; i++) {
                if (this.toString().equals(choices[i].toString())) {
                    return choices[i];
                }
            }
            throw new RuntimeException("unknown PhaseEvent type: "+this);
        }

    }
}