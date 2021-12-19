/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;

/**
 * Class that holds a set of atom pairs.  Takes a list of atoms in its
 * constructor, and forms an instance of AtomPair for each pair formed in the
 * list.  Each AtomPair instance can be accessed via the getSPair method. 
 */
public class AtomPairSet implements java.io.Serializable {

    /**
     * Constructor for AtomPairSet.
     * @param list The list of atoms for which the set of pairs is formed.
     */
    public AtomPairSet(IMoleculeList list) {
        aPairs = new IMolecule[list.size()-1][][];
        setAtoms(list);
    }
    
    /**
     * Returns atom pair for ith and jth atoms in set.
     */
    public IMolecule[] getAPair(int i, int j) {
        if(i==j) throw new IllegalArgumentException("Error: asking for pair formed with both atoms the same");
        return i<j ? aPairs[i][j-i-1] : aPairs[j][i-j-1];
    }

    private void setAtoms(IMoleculeList list) {
        int N = list.size();
        for(int i=0; i<N-1; i++) {
            aPairs[i] = new IMolecule[N-1-i][];
            for(int j=0; j<N-1-i; j++) {
                aPairs[i][j] = new IMolecule[2];
                aPairs[i][j][0] = list.get(i);
                aPairs[i][j][1] = list.get(i+j+1);
                if(aPairs[i][j][0].getType().getIndex()>aPairs[i][j][1].getType().getIndex()){//only for mixtures
                    aPairs[i][j][0] = aPairs[i][j][1];
                    aPairs[i][j][1] = list.get(i);
                }
            }
        }
    }
    
    private static final long serialVersionUID = 1L;
    private final IMolecule[][][] aPairs;
}
