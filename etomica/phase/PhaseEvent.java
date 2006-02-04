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
            new Type("Atom added"),
            new Type("Atom removed"),
            new Type("Atom changed index"),
            new Type("Max global index decreased")};
        public final EnumeratedType[] choices() {return CHOICES;}
    }
    public static final Type ATOM_ADDED =        Type.CHOICES[0];
    public static final Type ATOM_REMOVED =      Type.CHOICES[1];
    public static final Type ATOM_CHANGE_INDEX = Type.CHOICES[2];
    public static final Type GLOBAL_INDEX =      Type.CHOICES[3];
}//end of PhaseEvent
    