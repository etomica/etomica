package etomica.lattice.crystal;
import etomica.space3d.Vector3D;

/**
 * An 8-atom basis that for a diamond crystal.  Diamond is 4 fcc sites each 
 * with 2 subsites.
 *
 * @author David Kofke
 */
public class BasisCubicDiamond extends Basis {
    
    /**
     * Makes a diamond-on-fcc 8-atom basis.
     */
    public BasisCubicDiamond() {
        super(scaledPositions);
    }
    
    private static final Vector3D[] scaledPositions = new Vector3D[] {
        new Vector3D(0.00, 0.00, 0.00),
        new Vector3D(0.00, 0.50, 0.50),
        new Vector3D(0.50, 0.50, 0.00),
        new Vector3D(0.50, 0.00, 0.50),
        new Vector3D(0.25, 0.25, 0.25),
        new Vector3D(0.25, 0.75, 0.75),
        new Vector3D(0.75, 0.75, 0.25),
        new Vector3D(0.75, 0.25, 0.75)
    };

    private static final long serialVersionUID = 1L;
}