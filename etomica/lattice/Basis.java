package etomica.lattice;
import etomica.*;

/**
 * Atom factory that makes a basis for use on a lattice.
 *
 * @author David Kofke
 */
 
 /* History
  * 09/16/02 (DAK) new
  */
 
public abstract class Basis extends AtomFactory {
    
    public Basis(Simulation sim) {
        super(sim);
    }
    
    //add methods for neighbor distances
    //and other things
    
}//end of Basis