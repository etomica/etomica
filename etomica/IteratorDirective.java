package etomica;

public class IteratorDirective /*extends Constants.TypedConstant*/ {
    
    private Atom atom;
    private Phase phase;
    
//    protected IteratorDirective(String label) {super(label);}

/*    public static final IteratorDirective[] CHOICES = new IteratorDirective[] {
        new IteratorDirective("UP"),
        new IteratorDirective("DOWN"),
    };
    
    public final Constants.TypedConstant[] choices() {return CHOICES;}
    
    public static final IteratorDirective UP = CHOICES[0];
    public static final IteratorDirective DOWN = CHOICES[1];
 */   
    public final IteratorDirective setAtom(Atom a) {
        atom = a;
        return this;
    }
    public final Atom getAtom() {return atom;}
    
    public final IteratorDirective setPhase(Phase p) {
        phase = p;
        return this;
    }
    public final Phase getPhase() {return phase;}
    
    public static class Up extends IteratorDirective {
//        public Up() {super("Up");}
    }
    
    public static class Down extends IteratorDirective {
//        public Down() {super("Down");}
    }
    
}//end of IteratorDirective    