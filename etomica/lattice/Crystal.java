package etomica.lattice;
import etomica.AtomFactory;

/**
 * A class packaging together a Primitive and a Basis.
 */
 
public class Crystal {
    
    /**
     * Constructs instance using a simple single-site basis.
     * @param primitive
     */
    public Crystal(Primitive primitive) {
        this(primitive, new Basis(primitive.space));
    }
    /**
     * Constructs instance using the given primitive and basis.
     * @param primitive
     * @param basis
     */
    public Crystal(Primitive primitive, Basis basis) {
        this.primitive = primitive;
        this.basis = basis;
    }
    
    /**
     * Indicates crystal with one site at each primitive-vector combination,
     * with Atom at each site made by the given factory.
     * @param primitive
     * @param factory
     */
    public Crystal(Primitive primitive, AtomFactory factory) {
    	this(primitive, new Basis(primitive.space, 1, factory));
    }
    
    public Primitive getPrimitive() {return primitive;}
    
    public Basis getBasis() {return basis;}
    
    public String toString() {return primitive.toString();}
    
    protected Primitive primitive;
    protected Basis basis;
}//end of Crystal