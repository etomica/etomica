/*
 * Created on Mar 2, 2005
 */
package etomica.nbr;

import etomica.AtomPair;

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
    public boolean accept(AtomPair pair) {
        int diff = pair.atom0.node.getOrdinal() - pair.atom1.node.getOrdinal();
        if (isBonded != (diff == 1 || diff == -1) 
                || (!pair.atom0.inSameMolecule(pair.atom1))) {
            return false;
        }
        return subCriterion.accept(pair);
    }
    
    private boolean isBonded;

}
