/*
 * Created on Mar 2, 2005
 *
 * TODO To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package etomica.nbratom;

import etomica.AtomPair;
import etomica.nbr.NeighborCriterion;

/**
 * @author andrew
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
public class CriterionSpecies extends CriterionAdapter {

    public CriterionSpecies(NeighborCriterion criterion) {
        super(criterion);
    }
    
    public void setIntraSpecies(boolean b) {
        isIntraSpecies = b;
    }
    public boolean isIntraSpecies() {
        return isIntraSpecies;
    }
    
    public boolean accept(AtomPair pair) {
        if (isIntraSpecies != (pair.atom0.inSameSpecies(pair.atom1))) {
            return false;
        }
        return subCriterion.accept(pair);
    }
    
    private boolean isIntraSpecies;

}
