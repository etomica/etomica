package etomica.lattice;

import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.Primitive;
import etomica.space.Space;
import etomica.space.Vector;


/**
 * A lattice with sites given by the "atom" sites of a crystal.  Sites of
 * this lattice are instances of Space.Vector.  The dimension of a 
 * BravaisLatticeCrystal is one more than the dimension of the underlying Bravais
 * lattice forming the crystal; the extra index specifies the basis atom at
 * the site referenced by the other indices.
 */
public class BravaisLatticeCrystal extends BravaisLattice {

    /**
     * Constructs a lattice having sites given by the "atom" sites
     * of the given crystal.
     */
    public BravaisLatticeCrystal(Primitive primitive, Basis basis) {
        super(primitive);
        this.basis= basis;
        crystalIndex = new int[primitive.D];
        D = primitive.D + 1;
        position = getSpace().makeVector();
    }

    public Space getSpace() {
        return primitive.space;
    }

    /**
     * Returns the spatial dimension + 1.  The extra index specifies
     * the basis atom.
     */
    public int D() {
        return D;
    }

    /**
     * Returns a Vector instance giving the location of the referenced
     * site.  The first D-1 indices indicate the Bravais-lattice position, and 
     * the last index specifies the basis atom at the Bravais site.  The same
     * Vector instance is returned with each call.
     */
    public Object site(int[] index) {
        if(index.length != D) throw new IllegalArgumentException("index given to site method of lattice must have number of elements equal to dimension of lattice");
        System.arraycopy(index, 0, crystalIndex, 0, D-1);
        Vector latticePosition = (Vector)super.site(crystalIndex);
        Vector basisCoordinate = basis.getScaledCoordinates()[index[D-1]];
        Vector[] primitiveVectors = primitive.vectors();
        position.E(latticePosition);
        for (int i=0; i<basisCoordinate.D(); i++) {
            position.PEa1Tv1(basisCoordinate.x(i),primitiveVectors[i]);
        }
        return position;
    }
    
    public Basis getBasis() {
        return basis;
    }

    private static final long serialVersionUID = 1L;
    protected final Basis basis;
    private final int[] crystalIndex;
    private final int D;
    private final Vector position;
}
