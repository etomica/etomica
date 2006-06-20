package etomica.nbr;

import etomica.atom.AtomPair;
import etomica.atom.AtomSet;

/**
 * Pair criterion that judges whether two atoms are or are not in
 * the same molecule. Configurable to accept intra- or inter-molecular
 * pairs. 
 */

/*
 * Created on Mar 2, 2005
 */
public class CriterionInterMolecular extends CriterionAdapter {

    /**
     * Constructs criterion in default state such that intramolecular pairs are rejected,
     * and intermolecular pairs are accepted.
     */
    public CriterionInterMolecular(NeighborCriterion criterion) {
        super(criterion);
    }
   
    /**
     * Configures to accept intra- (if argument is true) or inter-
     * (if argument is false) molecular pairs.  Default is false.
     */
    public void setIntraMolecularCriterion(NeighborCriterion criterion) {
        intraCriterion = criterion;
    }

    /**
     * Flag indicating whether to accept intra- (if argument is true) or inter-
     * (if argument is false) molecular pairs.
     */
    public NeighborCriterion getIntraMolecularCriterion() {
        return intraCriterion;
    }
    
    /**
     * Returns false if pair is/isn't in same molecule (depending on setting
     * of intraMolecular); if matches this criterion, return value will be
     * that given by any subCriterion.
     */
    public boolean accept(AtomSet pair) {
        if ((((AtomPair)pair).atom0.inSameMolecule(((AtomPair)pair).atom1)) &&
                (intraCriterion == null || !intraCriterion.accept(pair))) {
            return false;
        }
        return subCriterion.accept(pair);
    }
    
    private NeighborCriterion intraCriterion;
}
