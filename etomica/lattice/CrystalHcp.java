package etomica.lattice;
import etomica.Simulation;
import etomica.Default;

/**
 * Hexagonal primitive with a 2-site hcp basis.
 */

 /* History
  * 09/27/02 (DAK) new
  */
public class CrystalHcp extends Crystal {
    
    public CrystalHcp(Simulation sim) {
        super(new PrimitiveHexagonal(sim));
        ((PrimitiveHexagonal)primitive).setA(Default.ATOM_SIZE);
        ((PrimitiveHexagonal)primitive).setC(Math.sqrt(8.0/3.0)*Default.ATOM_SIZE);
        siteFactory = new BasisHcp(sim, (PrimitiveHexagonal)primitive);
    }
    
    public String toString() {return "Hcp";}
    
}//end of CrystalHcp