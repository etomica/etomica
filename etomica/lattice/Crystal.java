package etomica.lattice;

import etomica.*;

/**
 * A class packaging together a Primitive and a Basis.
 */
 
public class Crystal {
    
    public Crystal(Primitive primitive) {
        this(primitive, new AtomFactoryMono(primitive.parentSimulation()));
    }
    public Crystal(Primitive primitive, AtomFactory factory) {
        this.primitive = primitive;
        this.siteFactory = factory;
    }
    
    public Primitive getPrimitive() {return primitive;}
    
    public AtomFactory getSiteFactory() {return siteFactory;}
    
    public String toString() {return primitive.toString();}
    
    protected Primitive primitive;
    protected AtomFactory siteFactory;
}//end of Crystal