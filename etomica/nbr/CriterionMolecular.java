package etomica.nbr;

import etomica.AtomPair;

/**
 * Pair criterion that judges whether two atoms are or are not in
 * the same molecule. Configurable to accept intra- or inter-molecular
 * pairs. 
 */

/*
 * Created on Mar 2, 2005
 */
public class CriterionMolecular extends CriterionAdapter {

    /**
     * Constructs criterion in default state such that intramolecular pairs are rejected,
     * and intermolecular pairs are accepted.
     */
    public CriterionMolecular(NeighborCriterion criterion) {
        super(criterion);
    }
   
    /**
     * Configures to accept intra- (if argument is true) or inter-
     * (if argument is false) molecular pairs.  Default is false.
     */
    public void setIntraMolecular(boolean b) {
        isIntraMolecular = b;
    }

    /**
     * Flag indicating whether to accept intra- (if argument is true) or inter-
     * (if argument is false) molecular pairs.
     */
    public boolean isIntraMolecular() {
        return isIntraMolecular;
    }
    
    /**
     * Returns false if pair is/isn't in same molecule (depending on setting
     * of intraMolecular); if matches this criterion, return value will be
     * that given by any subCriterion.
     */
    public boolean accept(AtomPair pair) {
        if (isIntraMolecular != (pair.atom0.inSameMolecule(pair.atom1))) {
            return false;
        }
        return subCriterion.accept(pair);
    }
    
    private boolean isIntraMolecular;

}
