package etomica.lattice;
import etomica.Space;
import etomica.Default;
import etomica.AtomFactory;

/**
 * Cubic primitive with a 4-site fcc basis.
 * 
 */

 /* History
  * 12/09/02 (skkwak) built another contructor to pass an instance of AtomFactory class into BasisCubicFcc 
  *                   to change atoms to sites (concrete sites under basis from BravisLattice) 
  * 09/23/02 (DAK) new
  */
public class CrystalFcc extends Crystal {

    public CrystalFcc(Space space) {
        super(new PrimitiveCubic(space));

        //set primitive to size for close packing if atoms are default size
        ((PrimitiveCubic)primitive).setSize(Math.sqrt(2.0)*Default.ATOM_SIZE);
        
        siteFactory = new BasisCubicFcc(space, (PrimitiveCubic)primitive);
    }

    public CrystalFcc(Space space, AtomFactory factory) {
        super(new PrimitiveCubic(space), factory);
        //set primitive to size for close packing if atoms are default size
        ((PrimitiveCubic)primitive).setSize(Math.sqrt(2.0)*Default.ATOM_SIZE);
        
        siteFactory = new BasisCubicFcc(space, factory, (PrimitiveCubic)primitive);
    }    
    
    public String toString() {return "Fcc";}
    
}//end of CrystalFcc