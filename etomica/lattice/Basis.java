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
 
 //maybe set up as an interface?
 
 // ? set a constructor in BravaisLattice to take a basis
 //   as an argument; add a getFactory method to the Basis
 //   interface; some implementations would extend AtomFactory
 //   and return 'this';
public abstract class Basis extends AtomFactory {
    
    public Basis(Space space, AtomSequencer.Factory seqFactory) {
        super(space, seqFactory);
    }
    
    //add methods for neighbor distances
    //and other things
    
}//end of Basis