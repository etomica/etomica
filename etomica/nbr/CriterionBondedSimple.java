/*
 * Created on Mar 2, 2005
 */
package etomica.nbr;

import etomica.atom.AtomPair;
import etomica.atom.AtomSet;

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
    public boolean accept(AtomSet pair) {
        int diff = ((AtomPair)pair).atom0.node.getIndex() - ((AtomPair)pair).atom1.node.getIndex();
        if (isBonded != (diff == 1 || diff == -1) 
                || (!((AtomPair)pair).atom0.inSameMolecule(((AtomPair)pair).atom1))) {
            return false;
        }
        return subCriterion.accept(pair);
    }
    
    private boolean isBonded;

}
