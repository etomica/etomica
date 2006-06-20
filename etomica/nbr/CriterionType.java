
package etomica.nbr;

import etomica.atom.Atom;
import etomica.atom.AtomSet;
import etomica.atom.AtomType;

/**
 * Filters atoms pairs to match a given pair of species.
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
     * Returns true if the species for the pair of atoms match 
     * the species given at construction (without regard to the
     * order of the pair), and if the wrapped criterion accept
     * also returns true.
     */
    public boolean accept(AtomSet atom) {
        AtomType atomType = ((Atom)atom).type;
        if (atomType == type) {
            return subCriterion.accept(atom);
        }
        return false;
    }
    
    public AtomType[] getTypes() {
        return new AtomType[]{type};
    }
    
    private final AtomType type;
}
