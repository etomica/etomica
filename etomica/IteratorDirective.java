package etomica;

/**
 * Encapsulation of a set of instructions that an AtomIterator
 * uses to select the atoms it presents on iteration.
 *
 * @author David Kofke
 */
public class IteratorDirective implements java.io.Serializable {
    
    private Atom atom1, atom2;
    private Direction direction;
    private int atomCount;
//    private Bounds bounds;
    
    public IteratorDirective() {
        this(UP);
    }
    public IteratorDirective(Direction direction) {
        set(direction);
        set();
//        set(ALL);
    }
    
    //returns itself as a convenience, so that it may be set while being passed as an
    //argument to a method
    public final IteratorDirective set() {
        atom1 = atom2 = null;
        atomCount = 0;
        return this;
    }
    public final IteratorDirective set(Atom a) {
        atom1 = a;
        atom2 = null;
        atomCount = 1;
//        bounds = FIRST;
        return this;
    }
    public final IteratorDirective set(Atom a1, Atom a2) {
        atom1 = a1;
        atom2 = a2;
        atomCount = 2;
//        bounds = FIRST_LAST;
        return this;
    }
    public final IteratorDirective set(Direction direction) {
        this.direction = direction;
        return this;
    }
/*    public final IteratorDirective set(Bounds bounds) {
        this.bounds = bounds;
        return this;
    }
*/    
//    public final Bounds bounds() {return bounds;}
    public final int atomCount() {return atomCount;}
    public final Direction direction() {return direction;}
    
    public final Atom atom1() {return atom1;}
    public final Atom atom2() {return atom2;}
    
    //IteratorDirective.Direction
    public static final class Direction extends Constants.TypedConstant {
            
        private Direction(String label) {super(label);}
        public static final Direction[] CHOICES = new Direction[] {
            new Direction("Up"),
            new Direction("Down"),
            new Direction("Singlet"),
        };
        
        public final Constants.TypedConstant[] choices() {return CHOICES;}
    }//end of Direction
    public static final Direction UP = Direction.CHOICES[0];
    public static final Direction DOWN = Direction.CHOICES[1];
    public static final Direction SINGLET = Direction.CHOICES[2];
    
    
    //IteratorDirective.Bounds
/*    public static final class Bounds extends Constants.TypedConstant {
            
        private Bounds(String label) {super(label);}
        public static final Bounds[] CHOICES = new Bounds[] {
            new Bounds("First"),
            new Bounds("FirstLast"),
            new Bounds("All")
        };
        
        public final Constants.TypedConstant[] choices() {return CHOICES;}
    }//end of Bounds
    public static final Bounds FIRST = Bounds.CHOICES[0];
    public static final Bounds FIRST_LAST = Bounds.CHOICES[1];
    public static final Bounds ALL = Bounds.CHOICES[2];
 */   
}//end of IteratorDirective    