package etomica.lattice.crystal;
import etomica.space3d.Vector3D;

/**
 * A 2-atom basis for a bcc crystal
 *
 * @author David Kofke
 */
 public class BasisCubicBcc extends Basis {
    
    /**
     * Makes a bcc 2-atom basis.
     */
    public BasisCubicBcc() {
        super(scaledPositions);
    }
    
    private static final Vector3D[] scaledPositions = new Vector3D[] {
            new Vector3D(0.0, 0.0, 0.0),
            new Vector3D(0.5, 0.5, 0.5)
    };
    
    private static final long serialVersionUID = 1L;
}