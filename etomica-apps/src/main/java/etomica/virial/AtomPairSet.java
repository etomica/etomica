/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.atom.IMoleculeList;
import etomica.atom.MoleculePair;

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
        aPairs = new MoleculePair[list.getMoleculeCount()-1][];
        setAtoms(list);
    }
    
    /**
     * Returns atom pair for ith and jth atoms in set.
     */
    public IMoleculeList getAPair(int i, int j) {
        if(i==j) throw new IllegalArgumentException("Error: asking for pair formed with both atoms the same");
        return i<j ? aPairs[i][j-i-1] : aPairs[j][i-j-1];
    }

    private void setAtoms(IMoleculeList list) {
        int N = list.getMoleculeCount();
        for(int i=0; i<N-1; i++) {
            aPairs[i] = new MoleculePair[N-1-i];
            for(int j=0; j<N-1-i; j++) {
                MoleculePair aPair = new MoleculePair();
                aPair.atom0 = list.getMolecule(i);
                aPair.atom1 = list.getMolecule(i+j+1);
                if(aPair.atom0.getType().getIndex()>aPair.atom1.getType().getIndex()){//only for mixtures
                	aPair.atom0 = aPair.atom1;
                	aPair.atom1 = list.getMolecule(i);
                }
                aPairs[i][j] = aPair;
            }
        }
    }
    
    private static final long serialVersionUID = 1L;
    private final MoleculePair[][] aPairs;
}
