package etomica.lattice;
import etomica.util.EnumeratedType;

public class LatticeEvent extends java.util.EventObject {
    
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
    
    public static final Type REBUILD = new Type("Rebuild");
    public static final Type ALL_SITE = new Type("All-site change");
    public static final Type RESET_NBRS = new Type("Reset neighbors");

    public static class Type extends EnumeratedType {
        protected Type(String label) {super(label);}
        public static Type[] choices() {
            return new Type[] {REBUILD,ALL_SITE,RESET_NBRS};
        }
    }
    
}//end of LatticeEvent