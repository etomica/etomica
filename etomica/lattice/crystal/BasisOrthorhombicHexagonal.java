package etomica.lattice.crystal;
import etomica.space2d.Vector2D;

/**
 * A 4-atom basis fot an fcc crystal.
 *
 * @author David Kofke
 */
public class BasisOrthorhombicHexagonal extends Basis {
    
    /**
     * Makes a fcc 4-atom basis.
     */
    public BasisOrthorhombicHexagonal() {
        super(scaledPositions);
    }
    
    private static final Vector2D[] scaledPositions = new Vector2D[] {
            new Vector2D(0.0, 0.0),
            new Vector2D(0.5, 0.5)
    };
    
    private static final long serialVersionUID = 1L;
}