package etomica.lattice;
import etomica.Space;
import etomica.Default;

/**
 * Cubic primitive with a 4-site fcc basis, on which each site 
 * is a 2-site diamond basis.
 */

 /* History
  * 09/26/02 (DAK) new
  */
public class CrystalDiamond extends Crystal {
    
    public CrystalDiamond(Space space) {
        super(new PrimitiveCubic(space));

        //set primitive to size for densest packing if atoms are default size
        ((PrimitiveCubic)primitive).setSize(4.0/Math.sqrt(3.0)*Default.ATOM_SIZE);
        
        //factory that makes the 2-atom sub-basis
        etomica.AtomFactory subsiteFactory = new BasisCubicFccDiamond(space, (PrimitiveCubic)primitive);

        //factory that makes the 4-atom fcc-on-cubic basis
        siteFactory = new BasisCubicFcc(space, subsiteFactory, (PrimitiveCubic)primitive);
    }
    
    public String toString() {return "Diamond";}
    
}//end of CrystalFcc