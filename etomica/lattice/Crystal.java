package etomica.lattice;

/**
 * A class packaging together a Primitive and a Basis.
 */
 
public class Crystal implements AbstractLattice {
    
    public Crystal(BravaisLattice lattice, Basis basis) {
        this.basis = basis;
        this.lattice = lattice;
    }
    
    public int D() {
        return lattice.D() + 1;
    }
    
    /* (non-Javadoc)
     * @see etomica.lattice.AbstractLattice#getSize()
     */
    public int[] getSize() {
        // TODO Auto-generated method stub
        return null;
    }
    /* (non-Javadoc)
     * @see etomica.lattice.AbstractLattice#site(int[])
     */
    public Object site(int[] index) {
        // TODO Auto-generated method stub
        return null;
    }
    /* (non-Javadoc)
     * @see etomica.lattice.AbstractLattice#sites()
     */
    public Object[] sites() {
        // TODO Auto-generated method stub
        return null;
    }
    public Basis getBasis() {return basis;}
    
    public BravaisLattice getLattice() {return lattice;}

    protected Basis basis;
    protected BravaisLattice lattice;
}//end of Crystal