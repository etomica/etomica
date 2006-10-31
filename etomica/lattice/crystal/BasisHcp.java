package etomica.lattice.crystal;
import etomica.space3d.Vector3D;

/**
 * A 2-atom basis for an hcp crystal.
 *
 * @author David Kofke
 */
public class BasisHcp extends Basis {
    
    public BasisHcp() {
        super(scaledPositions);
    }
    
    private static final Vector3D[] scaledPositions = new Vector3D[] {
        new Vector3D(0.0, 0.0, 0.0),
        new Vector3D(1.0/3.0, 1.0/3.0, 0.5),
    };

    private static final long serialVersionUID = 1L;
}
