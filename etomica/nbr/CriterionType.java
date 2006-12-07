
package etomica.nbr;

import etomica.atom.Atom;
import etomica.atom.AtomSet;
import etomica.atom.AtomType;

/**
 * Filters atoms to match a given AtomType.
 * 
 * @author Andrew Schultz
 */
public class CriterionType extends CriterionAdapter {

    public CriterionType(NeighborCriterion criterion, 
            AtomType type) {
        super(criterion);
        this.type = type;
    }
    
    /**
     * Returns true if the AtomType of the atom matches the AtomType given at 
     * construction and if the wrapped criterion accept also returns true.
     */
    public boolean accept(AtomSet atom) {
        if (((Atom)atom).getType() == type) {
            return subCriterion.accept(atom);
        }
        return false;
    }
    
    /**
     * Returns the AtomType accepted by this criterion.
     */
    public AtomType getType() {
        return type;
    }
    
    private static final long serialVersionUID = 1L;
    private final AtomType type;
}
