package etomica.phase;

import etomica.atom.Atom;
import etomica.graphics.DisplayPhaseListener;
import etomica.util.EnumeratedType;

/**
 * Event that conveys some happening with respect to a phase or the things it contains.
 *
 * @see PhaseListener
 * @see DisplayPhaseListener
 */
public class PhaseEvent extends java.util.EventObject {
    
    public PhaseEvent(Object source) {
        this(source, null);
    }
    public PhaseEvent(Object source, Type t) {
        super(source);
        type = t;
    }
    
    public void setType(Type t) {type = t;}
    public Type type() {return type;}
    
    public void setPhase(Phase p) {
        phase = p;
    }
    public Phase getPhase() {
        return phase;
    }
    
    public final PhaseEvent setAtom(Atom a) {atom = a; return this;}
    public Atom atom() {return atom;}
    
    public final void setIndex(int i) {index = i;}
    public final int getIndex() {return index;}
    
    protected Phase phase;
    protected Atom atom;
    protected Type type;
    protected int index;
    
    public static final Type ATOM_ADDED =        new Type("Atom added");
    public static final Type ATOM_REMOVED =      new Type("Atom removed");
    public static final Type ATOM_CHANGE_INDEX = new Type("Atom changed index");
    public static final Type GLOBAL_INDEX =      new Type("Max global index decreased");
    public static final Type PHASE_INFLATE =     new Type("Phase inflate");

    public static class Type extends EnumeratedType {

        protected Type(String label) {super(label);}

        public static Type[] choices() {
            return new Type[] {ATOM_ADDED,ATOM_REMOVED,ATOM_CHANGE_INDEX,GLOBAL_INDEX,PHASE_INFLATE};
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
}//end of PhaseEvent
    