package etomica.lattice.crystal;
import etomica.space3d.Vector3D;

/**
 * A 2-atom basis that makes a diamond crystal using a BravaisLattice
 * having a Cubic primitive with an fcc basis.  Each of the 4 fcc sites
 * is populated by a 2-atom basis arranged to yield the diamond structure.
 *
 * @author David Kofke
 */
 
 /* History
  * 09/26/02 (DAK) new
  * 01/19/04 (DAK) revised to extend Basis instead of AtomFactory
  */
 
public class BasisCubicFccDiamond extends BasisCubic {
    
    /**
     * Makes a diamond-on-fcc 2-atom basis using the given factory to make the atoms.
     */
    public BasisCubicFccDiamond(PrimitiveCubic primitive) {
        super(primitive, scaledPositions);
    }
    
    private static final Vector3D[] scaledPositions = new Vector3D[] {
            new Vector3D(0.0, 0.0, 0.0),
            new Vector3D(0.25, 0.25, 0.25),
    };
    
}//end of BasisCubicFcc