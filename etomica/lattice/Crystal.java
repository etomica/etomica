package etomica.lattice;

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
    
    public Primitive getPrimitive() {return primitive;}
    
    public Basis getBasis() {return basis;}
    
    public String toString() {return primitive.toString();}
    
    protected Primitive primitive;
    protected Basis basis;
}//end of Crystal