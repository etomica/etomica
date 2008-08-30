
package etomica.nbr;

import etomica.api.IAtomLeaf;
import etomica.api.IAtomSet;
import etomica.api.IAtomType;
import etomica.api.IAtomTypeLeaf;

/**
 * Filters atoms pairs to match a given pair of AtomTypes.
 * 
 * @author Andrew Schultz
 */
public class CriterionTypePair extends CriterionAdapter {

    public CriterionTypePair(NeighborCriterion criterion, 
            IAtomTypeLeaf type0, IAtomTypeLeaf type1) {
        super(criterion);
        this.type0 = type0;
        this.type1 = type1;
    }
    
    /**
     * Returns true if the AtomTypes for the pair of atoms match the AtomTypes 
     * given at construction (without regard to the order of the pair), and if 
     * the wrapped criterion also accepts the pair.
     */
    public boolean accept(IAtomSet pair) {
        IAtomType atom0Type = ((IAtomLeaf)pair.getAtom(0)).getType();
        IAtomType atom1Type = ((IAtomLeaf)pair.getAtom(1)).getType();
        if ( (atom0Type == type0 && atom1Type == type1) ||
             (atom0Type == type1 && atom1Type == type0) ) {
            return subCriterion.accept(pair);
        }
        return false;
    }
    
    /**
     * Returns the AtomTypes accepted by this NeighborCriterion
     */
    public IAtomType[] getTypes() {
        return new IAtomType[]{type0,type1};
    }
    
    private static final long serialVersionUID = 1L;
    private final IAtomTypeLeaf type0;
    private final IAtomTypeLeaf type1;
}
