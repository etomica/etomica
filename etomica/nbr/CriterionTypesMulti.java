
package etomica.nbr;

import etomica.api.IAtomSet;
import etomica.api.IAtomType;

/**
 * Filters AtomSets to match a given set of AtomTypes.  CriterionType and
 * CriterionTypePair should be used for single and pairs of Atoms.
 * 
 * @author Andrew Schultz
 */
public class CriterionTypesMulti extends CriterionAdapter {

    public CriterionTypesMulti(NeighborCriterion criterion, 
            IAtomType[] types) {
        super(criterion);
        this.types = (IAtomType[])types.clone();
    }
    
    /**
     * Returns true if the AtomTypes for the pair of atoms match the AtomTypes 
     * given at construction (without regard to the order of the AtomSet), and 
     * if the wrapped criterion also accepts the AtomSet.
     */
    public boolean accept(IAtomSet atoms) {
        final int nAtoms = atoms.getAtomCount();
        for (int i=0; i<types.length; i++) {
            boolean accepted = false;
            for (int j=0; j<nAtoms; j++) {
                if (atoms.getAtom(j).getType() == types[i]) {
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
    public IAtomType[] getTypes() {
        return (IAtomType[])types.clone();
    }
    
    private static final long serialVersionUID = 1L;
    private final IAtomType[] types;
}
