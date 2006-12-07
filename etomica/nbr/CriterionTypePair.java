
package etomica.nbr;

import etomica.atom.AtomPair;
import etomica.atom.AtomSet;
import etomica.atom.AtomType;

/**
 * Filters atoms pairs to match a given pair of AtomTypes.
 * 
 * @author Andrew Schultz
 */
public class CriterionTypePair extends CriterionAdapter {

    public CriterionTypePair(NeighborCriterion criterion, 
            AtomType type0, AtomType type1) {
        super(criterion);
        this.type0 = type0;
        this.type1 = type1;
    }
    
    /**
     * Returns true if the AtomTypes for the pair of atoms match the AtomTypes 
     * given at construction (without regard to the order of the pair), and if 
     * the wrapped criterion also accepts the pair.
     */
    public boolean accept(AtomSet pair) {
        AtomType atom0Type = ((AtomPair)pair).atom0.getType();
        AtomType atom1Type = ((AtomPair)pair).atom1.getType();
        if ( (atom0Type == type0 && atom1Type == type1) ||
             (atom0Type == type1 && atom1Type == type0) ) {
            return subCriterion.accept(pair);
        }
        return false;
    }
    
    /**
     * Returns the AtomTypes accepted by this NeighborCriterion
     */
    public AtomType[] getTypes() {
        return new AtomType[]{type0,type1};
    }
    
    private static final long serialVersionUID = 1L;
    private final AtomType type0;
    private final AtomType type1;
}
