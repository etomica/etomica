package etomica.lattice.crystal;
import etomica.Space3D;

/**
 * A 2-atom basis that makes a bcc crystal on a BravaisLattice
 * having a Cubic primitive.
 *
 * @author David Kofke
 */
 
public class BasisCubicBcc extends BasisCubic {
    
    /**
     * Makes a bcc 2-atom basis.
     * @param primitive Primitive of the cubic lattice housing this basis.
     * Needed to ensure that separation of basis atoms is consistent with
     * spacing of atoms on lattice.
     */
    public BasisCubicBcc(PrimitiveCubic primitive) {
        super(primitive, scaledPositions);
    }
    
    private static final Space3D.Vector[] scaledPositions = new Space3D.Vector[] {
            new Space3D.Vector(0.0, 0.0, 0.0),
            new Space3D.Vector(0.5, 0.5, 0.5)
    };
    
}//end of BasisCubicBcc