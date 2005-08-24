package etomica.lattice;
import etomica.SimulationEvent;
import etomica.util.Constants;
import etomica.util.EnumeratedType;

public class LatticeEvent extends SimulationEvent {
    
    protected AbstractLattice lattice;
    protected Type type;
    
    public LatticeEvent(AbstractLattice source) {
        this(source, null);
    }
    public LatticeEvent(AbstractLattice source, Type t) {
        super(source);
        lattice = source;
        type = t;
    }
    
    public void setType(Type t) {type = t;}
    public Type type() {return type;}
    
    public final LatticeEvent setLattice(AbstractLattice c) {lattice = c; return this;}
    public final AbstractLattice lattice() {return lattice;}
    
    public static class Type extends EnumeratedType {
        private Type(String label) {super(label);}
        public static final Type[] CHOICES = new Type[] {
            new Type("Rebuild"),
            new Type("All-site change"),
            new Type("Reset neighbors")};
        public final EnumeratedType[] choices() {return CHOICES;}
    }
    public static final Type REBUILD = Type.CHOICES[0];
    public static final Type ALL_SITE = Type.CHOICES[1];
    public static final Type RESET_NBRS = Type.CHOICES[2];
    
}//end of LatticeEvent