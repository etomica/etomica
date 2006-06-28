package etomica.nbr;

import etomica.atom.AtomPair;
import etomica.atom.AtomSet;

/**
 * Pair criterion that judges whether two atoms are or are not in the same 
 * molecule.  Intermolecular pairs are always accepted.  An optional 
 * intramolecular may be used to filter out intramolecular pairs.
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
     * Configures to use the given criterion for intramolecular pairs.
     */
    public void setIntraMolecularCriterion(NeighborCriterion criterion) {
        intraCriterion = criterion;
    }

    /**
     * Returns the intramolecular criterion, or null if none is in use.
     */
    public NeighborCriterion getIntraMolecularCriterion() {
        return intraCriterion;
    }
    
    public boolean accept(AtomSet pair) {
        // Only ask the intracriterion if it exists and the pair is intramolecular. 
        if (intraCriterion != null && (((AtomPair)pair).atom0.inSameMolecule(((AtomPair)pair).atom1)) &&
                !intraCriterion.accept(pair)) {
            return false;
        }
        return subCriterion.accept(pair);
    }
    
    private NeighborCriterion intraCriterion;
}
