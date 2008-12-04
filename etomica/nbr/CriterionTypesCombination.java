
package etomica.nbr;

import etomica.api.IAtomLeaf;
import etomica.api.IAtomList;
import etomica.api.IAtomType;
import etomica.api.IAtomTypeLeaf;

/**
 * Filters AtomSets to match a given set of AtomTypes.  CriterionType and
 * CriterionTypePair should be used for single and pairs of Atoms.
 * 
 * @author Andrew Schultz
 */
public class CriterionTypesCombination extends CriterionAdapter {

    public CriterionTypesCombination(NeighborCriterion criterion, 
            IAtomTypeLeaf[] types) {
        super(criterion);
        this.types = types.clone();
    }
    
    /**
     * Returns true if the AtomTypes for the pair of atoms match the AtomTypes 
     * given at construction (without regard to the order of the AtomSet), and 
     * if the wrapped criterion also accepts the AtomSet.
     */
    public boolean accept(IAtomList atoms) {
        final int nAtoms = atoms.getAtomCount();
        for (int i=0; i<nAtoms; i++) {
            boolean accepted = false;
            for (int j=0; j<types.length; j++) {
                if (((IAtomLeaf)atoms.getAtom(i)).getType() == types[j]) {
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
    public IAtomTypeLeaf[] getTypes() {
        return types;
    }
    
    private static final long serialVersionUID = 1L;
    private final IAtomTypeLeaf[] types;
}
