package etomica.lattice;
import etomica.Space;
import etomica.Default;

/**
 * Hexagonal primitive with a 2-site hcp basis.
 */

 /* History
  * 09/27/02 (DAK) new
  */
public class CrystalHcp extends Crystal {
    
    public CrystalHcp(Space space) {
        super(new PrimitiveHexagonal(space));
        ((PrimitiveHexagonal)primitive).setA(Default.ATOM_SIZE);
        ((PrimitiveHexagonal)primitive).setC(Math.sqrt(8.0/3.0)*Default.ATOM_SIZE);
        siteFactory = new BasisHcp(space, (PrimitiveHexagonal)primitive);
    }
    
    public String toString() {return "Hcp";}
    
}//end of CrystalHcp