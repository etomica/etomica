package etomica.lattice;
import etomica.Space;
import etomica.Default;

/**
 * Cubic primitive with a 2-site bcc basis.
 */

 /* History
  * 09/23/02 (DAK) new
  */
public class CrystalBcc extends Crystal {
    
    public CrystalBcc(Space space) {
        super(new PrimitiveCubic(space));
        ((PrimitiveCubic)primitive).setSize(2.0/Math.sqrt(3.0)*Default.ATOM_SIZE);
        siteFactory = new BasisCubicBcc(space, (PrimitiveCubic)primitive);
    }
    
    public String toString() {return "Bcc";}
    
}//end of CrystalBcc