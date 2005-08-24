package etomica.phase;

import etomica.SimulationEvent;
import etomica.atom.Atom;
import etomica.graphics.DisplayPhaseListener;
import etomica.space.Vector;
import etomica.util.Constants;
import etomica.util.EnumeratedType;

/**
 * Event that conveys some happening with respect to a phase or the things it contains.
 *
 * @see PhaseListener
 * @see DisplayPhaseListener
 */
public class PhaseEvent extends SimulationEvent {
    
    protected Phase phase;
    protected Atom atom;
    protected Vector point;
    protected Type type;
    
    //for boundary inflation event
    public double isoScale;
    public Vector anisoScale;
    public boolean isotropic;
    
    public PhaseEvent(Object source) {
        this(source, null);
    }
    public PhaseEvent(Object source, Type t) {
        super(source);
        type = t;
    }
    
    public void setType(Type t) {type = t;}
    public Type type() {return type;}
    
    public final PhaseEvent setPhase(Phase p) {phase = p; return this;}
    public final Phase phase() {return phase;}
    
    public final PhaseEvent setPoint(Vector p) {point = p; return this;}
    public Vector point() {return point;}
    
    public final PhaseEvent setAtom(Atom a) {atom = a; return this;}
    public Atom atom() {return atom;}
    
    public final PhaseEvent setScale(double s) {isoScale = s; isotropic = true; return this;}
    public final PhaseEvent setScale(Vector s) {anisoScale = s; isotropic = false; return this;}
    
    public static class Type extends EnumeratedType {
        private Type(String label) {super(label);}
        public static final Type[] CHOICES = new Type[] {
            new Type("Point selected"),
            new Type("Atom added"),
            new Type("Atom removed"),
            new Type("Atom selected"),
            new Type("Atom released"),
            new Type("Boundary PhaseInflate"),
            new Type("Reset")};
        public final EnumeratedType[] choices() {return CHOICES;}
    }
    public static final Type POINT_SELECTED =   Type.CHOICES[0];
    public static final Type ATOM_ADDED =       Type.CHOICES[1];
    public static final Type ATOM_REMOVED =     Type.CHOICES[2];
    public static final Type ATOM_SELECTED =    Type.CHOICES[3];
    public static final Type ATOM_RELEASED =    Type.CHOICES[4];
    public static final Type BOUNDARY_INFLATE = Type.CHOICES[5];
    public static final Type RESET =            Type.CHOICES[6];
}//end of PhaseEvent
    