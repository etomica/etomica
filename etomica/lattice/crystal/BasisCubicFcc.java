package etomica.lattice.crystal;
import etomica.space3d.Vector3D;

/**
 * A 4-atom basis that makes an fcc crystal on a BravaisLattice
 * having a Cubic primitive.
 *
 * @author David Kofke
 */
 
public class BasisCubicFcc extends BasisCubic {
    
    /**
     * Makes a fcc 4-atom basis using the given factory to make the atoms.
     * @param primitive Primitive of the cubic lattice housing this basis.
     * Needed to ensure that separation of basis atoms is consistent with
     * spacing of atoms on lattice.
     */
    public BasisCubicFcc(PrimitiveCubic primitive) {
        super(primitive, scaledPositions);
    }
    
    private static final Vector3D[] scaledPositions = new Vector3D[] {
            new Vector3D(0.0, 0.0, 0.0),
            new Vector3D(0.0, 0.5, 0.5),
            new Vector3D(0.5, 0.5, 0.0),
            new Vector3D(0.5, 0.0, 0.5)
    };
    
}//end of BasisCubicFcc