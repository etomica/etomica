package etomica.lattice;
import etomica.Simulation;

/**
 * Cubic primitive with a 2-site bcc basis.
 */

 /* History
  * 09/23/02 (DAK) new
  */
public class CrystalBcc extends Crystal {
    
    public CrystalBcc(Simulation sim) {
        super(new PrimitiveCubic(sim));
        siteFactory = new BasisCubicBcc(sim, (PrimitiveCubic)primitive);
    }
    
    public String toString() {return "Bcc (cubic)";}
    
}//end of CrystalBcc