package etomica.lattice;
import etomica.Space;
import etomica.Default;

/**
 * Cubic primitive with a 4-site fcc basis.
 */

 /* History
  * 09/23/02 (DAK) new
  */
public class CrystalFcc extends Crystal {
    
    public CrystalFcc(Space space) {
        super(new PrimitiveCubic(space));
        
        //set primitive to size for close packing if atoms are default size
        ((PrimitiveCubic)primitive).setSize(Math.sqrt(2.0)*Default.ATOM_SIZE);
        
        siteFactory = new BasisCubicFcc(space, (PrimitiveCubic)primitive);
    }
    
    public String toString() {return "Fcc";}
    
}//end of CrystalFcc