/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.IAtom;

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
