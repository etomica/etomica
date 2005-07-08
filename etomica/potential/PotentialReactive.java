package etomica.potential;

import etomica.Atom;

/**
 * Interface for potentials that change the bond state of an atom.
 */
public interface PotentialReactive {
        
    public BondChangeData[] getBondChangeData();
        
    public static class BondChangeData implements java.io.Serializable {
        public Atom atom;
        public Atom[] oldPartners;
        public Atom[] newPartners;
        public Atom getAtom() {return atom;}
        public Atom[] getOldPartners() {return oldPartners;}
        public Atom[] getNewPartners() {return newPartners;}
    }//end of BondChangeData
        
}//end of PotentialReactive
