package etomica.potential;

import etomica.api.IAtom;

/**
 * Interface for potentials that change the bond state of an atom.
 */
public interface PotentialReactive {
        
    public BondChangeData[] getBondChangeData();
        
    public static class BondChangeData implements java.io.Serializable {
        private static final long serialVersionUID = 1L;
        public IAtom atom;
        public IAtom[] oldPartners;
        public IAtom[] newPartners;
        public IAtom getAtom() {return atom;}
        public IAtom[] getOldPartners() {return oldPartners;}
        public IAtom[] getNewPartners() {return newPartners;}
    }//end of BondChangeData
        
}//end of PotentialReactive
