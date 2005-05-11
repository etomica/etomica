package etomica.lattice;

import etomica.Space;
import etomica.space.Vector;

/**
 * A class packaging together a lattice and a basis to define a crystal lattice.
 * Physically this would represent a molecular crystal, in which the basis describes
 * placement of the atoms of a molecule, with the molecule repeated on the regular site
 * of a lattice. The sites of this class are an array of vectors indicating the placement 
 * the "atoms" for the "molecule" at the Bravais lattice site for the given primitive.
 */
public class Crystal implements AbstractLattice {
    
    /**
     * Constructs a crystal using the given basis and a Bravais lattice 
     * based on the given primitive. 
     */
    public Crystal(Primitive primitive, Basis basis) {
        this(new BravaisLattice(primitive), basis);
    }
    
    public Crystal(BravaisLattice lattice, Basis basis) {
        this.basis = basis;
        this.lattice = lattice;
        this.space = lattice.getSpace();
        latticeIndex = new int[lattice.D()];
    }
    
    public int D() {
        return lattice.D();
    }
    
    public Space getSpace() {
        return space;
    }
    
    /**
     * Returns the basis at the Bravais-lattice site corresponding
     * to the given index.
     * 
     * @return a Space.Vector array with the positions of the basis sites
     */
    public Object site(int[] index) {
        Vector latticePosition = (Vector)lattice.site(index);
        Vector[] basisPositions = basis.positions();
        Vector[] positions = space.makeVectorArray(basis.size());
        for(int i=basis.size()-1; i>=0; i--) {
            positions[i].Ev1Pv2(latticePosition,basisPositions[i]);
        }
        return positions;
    }
    
    public Basis getBasis() {return basis;}
    
    public BravaisLattice getLattice() {return lattice;}

    protected final Basis basis;
    protected final BravaisLattice lattice;
    private final int[] latticeIndex;
    private final Space space;
}//end of Crystal