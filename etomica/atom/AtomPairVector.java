package etomica.atom;

import etomica.AtomPair;
import etomica.space.Vector;


/**
 * Atom pair that includes a vector giving the displacement 
 * applied to atom0 that makes it the nearest image of atom1.
 * Thus, nearest-image separation is 
 * atom1.coord.position - atom0.coord.position - nearestImageVector
 */

/*
 * History
 * Created on Feb 18, 2005 by kofke
 */
public class AtomPairVector extends AtomPair {

    public Vector nearestImageVector;
    
    
    /* (non-Javadoc)
     * @see etomica.AtomPair#copyTo(etomica.AtomPair)
     */
    public void copyTo(AtomPair pair) {
        super.copyTo(pair);
        ((AtomPairVector)pair).nearestImageVector = this.nearestImageVector;
    }
}
