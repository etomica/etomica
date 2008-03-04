package etomica.nbr;

import etomica.api.IAtomSet;
import etomica.atom.IAtomLeaf;

/**
 * @author andrew
 * Used for bonding potentials. The accept method returns true if the atoms are
 * adjacent in the list of atoms and they have the same parent. 
 */
public class CriterionBondedSimple extends CriterionAdapter {

    public CriterionBondedSimple(NeighborCriterion criterion) {
        super(criterion);
    }
    
    public void setBonded(boolean b) {
        isBonded = b;
    }
    public boolean isBonded() {
        return isBonded;
    }
    
    // always enforce intramolecularity
    public boolean accept(IAtomSet pair) {
        int diff = pair.getAtom(0).getIndex() - pair.getAtom(1).getIndex();
        if (isBonded != (diff == 1 || diff == -1) 
                || (((IAtomLeaf)pair.getAtom(0)).getParentGroup() != ((IAtomLeaf)pair.getAtom(1)).getParentGroup())) {
            return false;
        }
        return subCriterion.accept(pair);
    }
    
    private static final long serialVersionUID = 1L;
    private boolean isBonded;
}
