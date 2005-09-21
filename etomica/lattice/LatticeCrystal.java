package etomica.lattice;

import etomica.space.Space;
import etomica.space.Vector;


/**
 * A lattice with sites given by the "atom" sites of a crystal.  Sites of
 * this lattice are instances of Space.Vector (unlike the sites of a Crystal
 * which are an array of vectors, one for each atom in the basis).  The dimension
 * of a LatticeCrystal is one more than the dimension of the underlying Bravais
 * lattice forming the crystal; the extra index specifies the basis atom at
 * the site referenced by the other indices.
 */

/*
 * History
 * Created on Jan 5, 2005 by kofke
 */
public class LatticeCrystal implements SpaceLattice, java.io.Serializable {

    /**
     * Constructs a lattice having sites given by the "atom" sites
     * of the given crystal.
     */
    public LatticeCrystal(Crystal crystal) {
        this.crystal = crystal;
        crystalIndex = new int[crystal.D()];
        D = crystal.D() + 1;
        positions = getSpace().makeVector();
    }

    /* (non-Javadoc)
     * @see etomica.lattice.SpaceLattice#space()
     */
    public Space getSpace() {
        return crystal.getSpace();
    }

    /**
     * Equal to the spatial dimension + 1.  The extra index specifies
     * the basis atom.
     */
    public int D() {
        return D;
    }

    /**
     * Returns a etomica.space.Vector instance giving the location of the referenced
     * site.  The first D-1 indices indicate the Bravais-lattice position, and 
     * the last index specifies the basis atom at the Bravais site.  The same
     * Vector instance is returned with each call.
     */
    public Object site(int[] index) {
        if(index.length != D) throw new IllegalArgumentException("index given to site method of lattice must have number of elements equal to dimension of lattice");
        System.arraycopy(index, 0, crystalIndex, 0, D-1);
        Vector latticePosition = (Vector)crystal.getLattice().site(crystalIndex);
        Vector[] basisPositions = crystal.getBasis().positions();
        positions.Ev1Pv2(latticePosition,basisPositions[index[D-1]]);
        return positions;
    }
    
    public Crystal getCrystal() {
        return crystal;
    }
    
    protected final Crystal crystal;
    private final int[] crystalIndex;
    private final int D;
    private final Vector positions;
}
