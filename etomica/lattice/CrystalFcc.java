package etomica.lattice;
import etomica.Simulation;

/**
 * Cubic primitive with a 4-site fcc basis.
 */

 /* History
  * 09/23/02 (DAK) new
  */
public class CrystalFcc extends Crystal {
    
    public CrystalFcc(Simulation sim) {
        super(new PrimitiveCubic(sim));
        siteFactory = new BasisCubicFcc(sim, (PrimitiveCubic)primitive);
    }
    
    public String toString() {return "Fcc (cubic)";}
    
}//end of CrystalFcc