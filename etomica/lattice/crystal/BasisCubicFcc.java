package etomica.lattice.crystal;
import etomica.space3d.Vector3D;

/**
 * A 4-atom basis fot an fcc crystal.
 *
 * @author David Kofke
 */
public class BasisCubicFcc extends Basis {
    
    /**
     * Makes a fcc 4-atom basis.
     */
    public BasisCubicFcc() {
        super(scaledPositions);
    }
    
    private static final Vector3D[] scaledPositions = new Vector3D[] {
            new Vector3D(0.0, 0.0, 0.0),
            new Vector3D(0.0, 0.5, 0.5),
            new Vector3D(0.5, 0.5, 0.0),
            new Vector3D(0.5, 0.0, 0.5)
    };
    
    private static final long serialVersionUID = 1L;
}