package etomica.lattice;
import etomica.Simulation;
import etomica.Default;

/**
 * Cubic primitive with a 4-site fcc basis.
 */

 /* History
  * 09/23/02 (DAK) new
  */
public class CrystalFcc extends Crystal {
    
    public CrystalFcc(Simulation sim) {
        super(new PrimitiveCubic(sim));
        
        //set primitive to size for close packing if atoms are default size
        ((PrimitiveCubic)primitive).setSize(Math.sqrt(2.0)*Default.ATOM_SIZE);
        
        siteFactory = new BasisCubicFcc(sim, (PrimitiveCubic)primitive);
    }
    
    public String toString() {return "Fcc (cubic)";}
    
}//end of CrystalFcc