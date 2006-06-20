
package etomica.nbr;

import etomica.atom.AtomSet;
import etomica.atom.AtomType;

/**
 * Filters atoms pairs to match a given pair of species.
 * 
 * @author Andrew Schultz
 */
public class CriterionTypesMulti extends CriterionAdapter {

    public CriterionTypesMulti(NeighborCriterion criterion, 
            AtomType[] types) {
        super(criterion);
        this.types = (AtomType[])types.clone();
    }
    
    /**
     * Returns true if the species for the pair of atoms match 
     * the species given at construction (without regard to the
     * order of the pair), and if the wrapped criterion accept
     * also returns true.
     */
    public boolean accept(AtomSet atoms) {
        final int nAtoms = atoms.count();
        for (int i=0; i<types.length; i++) {
            boolean accepted = false;
            for (int j=0; j<nAtoms; j++) {
                if (atoms.getAtom(j).type == types[i]) {
                    accepted = true;
                }
            }
            if (!accepted) {
                return false;
            }
        }
        return subCriterion.accept(atoms);
    }
    
    public AtomType[] getTypes() {
        return (AtomType[])types.clone();
    }
    
    private final AtomType[] types;
}
