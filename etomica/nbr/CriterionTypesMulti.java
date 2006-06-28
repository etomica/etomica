
package etomica.nbr;

import etomica.atom.AtomSet;
import etomica.atom.AtomType;

/**
 * Filters AtomSets to match a given set of AtomTypes.  CriterionType and
 * CriterionTypePair should be used for single and pairs of Atoms.
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
     * Returns true if the AtomTypes for the pair of atoms match the AtomTypes 
     * given at construction (without regard to the order of the AtomSet), and 
     * if the wrapped criterion also accepts the AtomSet.
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
    
    /**
     * Returns the AtomTypes accepted by this NeighborCriterion
     */
    public AtomType[] getTypes() {
        return (AtomType[])types.clone();
    }
    
    private final AtomType[] types;
}
