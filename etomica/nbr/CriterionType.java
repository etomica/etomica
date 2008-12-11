
package etomica.nbr;

import etomica.api.IAtomList;
import etomica.api.IAtomTypeLeaf;

/**
 * Filters atoms to match a given AtomType.
 * 
 * @author Andrew Schultz
 */
public class CriterionType extends CriterionAdapter {

    public CriterionType(NeighborCriterion criterion, 
            IAtomTypeLeaf type) {
        super(criterion);
        this.type = type;
    }
    
    /**
     * Returns true if the AtomType of the atom matches the AtomType given at 
     * construction and if the wrapped criterion accept also returns true.
     */
    public boolean accept(IAtomList atom) {
        if (atom.getAtom(0).getType() == type) {
            return subCriterion.accept(atom);
        }
        return false;
    }
    
    /**
     * Returns the AtomType accepted by this criterion.
     */
    public IAtomTypeLeaf getType() {
        return type;
    }
    
    private static final long serialVersionUID = 1L;
    private final IAtomTypeLeaf type;
}
