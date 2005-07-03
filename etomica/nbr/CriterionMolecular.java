package etomica.nbr;

import etomica.AtomPair;

/**
 * @author andrew
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 */

/*
 * Created on Mar 2, 2005
 */
public class CriterionMolecular extends CriterionAdapter {

    public CriterionMolecular(NeighborCriterion criterion) {
        super(criterion);
    }
    
    public void setIntraMolecular(boolean b) {
        isIntraMolecular = b;
    }
    public boolean isIntraMolecular() {
        return isIntraMolecular;
    }
    
    public boolean accept(AtomPair pair) {
        if (isIntraMolecular != (pair.atom0.inSameMolecule(pair.atom1))) {
            return false;
        }
        return subCriterion.accept(pair);
    }
    
    private boolean isIntraMolecular;

}
