package etomica.lattice;

import etomica.Space;

/**
 * A class packaging together a Primitive and a Basis.
 */
public class Crystal implements AbstractLattice {
    
    public Crystal(SpaceLattice lattice, Basis basis) {
        this.basis = basis;
        this.lattice = lattice;
        this.space = lattice.space();
        latticeIndex = new int[lattice.D()];
    }
    
    public int D() {
        return lattice.D();
    }
    
    public Space space() {
        return space;
    }
    
    /**
     * Returns the basis at the Bravais-lattice site corresponding
     * to the given index.
     * 
     * @return a Space.Vector array with the positions of the basis sites
     */
    public Object site(int[] index) {
        Space.Vector latticePosition = (Space.Vector)lattice.site(index);
        Space.Vector[] basisPositions = basis.positions();
        Space.Vector[] positions = space.makeVectorArray(basis.size());
        for(int i=basis.size()-1; i>=0; i--) {
            positions[i].Ev1Pv2(latticePosition,basisPositions[i]);
        }
        return positions;
    }
    
    public Basis getBasis() {return basis;}
    
    public SpaceLattice getLattice() {return lattice;}

    protected final Basis basis;
    protected final SpaceLattice lattice;
    private final int[] latticeIndex;
    private final Space space;
}//end of Crystal