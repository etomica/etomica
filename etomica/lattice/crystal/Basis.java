package etomica.lattice.crystal;
import etomica.space.IVector;

/**
 * General basis class that hold scaled coordinates of atoms within a unit cell.
 *
 * @author David Kofke
 */
public class Basis implements java.io.Serializable {
    

    /**
     * @param scaledCoordinates basis coordinates for the case in which the
     * primitive lattice constant (getSize) is unity.  Given instance is kept
     * for interal representation of basis, so changes to scaledCoordinates
     * will affect the basis.
     */
    public Basis(IVector[] scaledCoordinates) {
        this.scaledCoordinates = scaledCoordinates;
    }
    
    /**
     * Returns scaled coordinates
     */
    public IVector[] getScaledCoordinates() {
        return scaledCoordinates;
    }
    
    private static final long serialVersionUID = 1L;
    private final IVector[] scaledCoordinates;
    
}
