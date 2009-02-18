package etomica.potential;

import etomica.api.IAtomLeaf;

/**
 * Interface for potentials that change the bond state of an atom.
 */
public interface PotentialReactive {
        
    public BondChangeData[] getBondChangeData();
        
    public static class BondChangeData implements java.io.Serializable {
        private static final long serialVersionUID = 1L;
        public IAtomLeaf atom;
        public IAtomLeaf[] oldPartners;
        public IAtomLeaf[] newPartners;
        public IAtomLeaf getAtom() {return atom;}
        public IAtomLeaf[] getOldPartners() {return oldPartners;}
        public IAtomLeaf[] getNewPartners() {return newPartners;}
    }//end of BondChangeData
        
}//end of PotentialReactive
