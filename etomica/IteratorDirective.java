package etomica;

public class IteratorDirective implements java.io.Serializable {
    
    private Atom firstAtom, lastAtom;
    private Phase phase;
    private Direction direction;
    private Bounds bounds;
    
    public IteratorDirective() {
        this(UP);
    }
    public IteratorDirective(Direction direction) {
        set(direction);
    }
    
    public final IteratorDirective set(Atom first) {
        firstAtom = first;
        bounds = FIRST;
        return this;
    }
    public final IteratorDirective set(Atom first, Atom last) {
        firstAtom = first;
        lastAtom = last;
        bounds = FIRST_LAST;
        return this;
    }
    public final IteratorDirective set(Direction direction) {
        this.direction = direction;
        return this;
    }
    public final IteratorDirective set(Bounds bounds) {
        this.bounds = bounds;
        return this;
    }
    
    public final Bounds bounds() {return bounds;}
    public final Direction direction() {return direction;}
    public final Atom firstAtom() {return firstAtom;}
    public final Atom lastAtom() {return lastAtom;}
    
    public final IteratorDirective setPhase(Phase p) {
        phase = p;
        return this;
    }
    public final Phase getPhase() {return phase;}
        
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
    
    //IteratorDirective.Direction
    public static final class Bounds extends Constants.TypedConstant {
            
        private Bounds(String label) {super(label);}
        public static final Bounds[] CHOICES = new Bounds[] {
            new Bounds("First"),
            new Bounds("FirstLast"),
            new Bounds("All")
        };
        
        public final Constants.TypedConstant[] choices() {return CHOICES;}
    }//end of Direction
    public static final Bounds FIRST = Bounds.CHOICES[0];
    public static final Bounds FIRST_LAST = Bounds.CHOICES[1];
    public static final Bounds ALL = Bounds.CHOICES[2];
    
}//end of IteratorDirective    