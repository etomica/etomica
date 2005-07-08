package etomica.virial;

import etomica.Atom;
import etomica.AtomPair;
import etomica.atom.AtomList;
import etomica.atom.iterator.AtomIteratorListTabbed;

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
    public AtomPairSet(AtomList list) {
        aPairs = new AtomPair[list.size()-1][];
        setAtoms(list);
    }
    
    /**
     * Returns atom pair for ith and jth atoms in set.
     */
    public AtomPair getAPair(int i, int j) {
        if(i==j) throw new IllegalArgumentException("Error: asking for pair formed with both atoms the same");
        return i<j ? aPairs[i][j-i-1] : aPairs[j][i-j-1];
    }

    private void setAtoms(AtomList list) {
        AtomIteratorListTabbed iterator = new AtomIteratorListTabbed(list);
        int N = list.size();
        Atom[] atoms = new Atom[N];
        iterator.reset();
        int k=0;
        while(iterator.hasNext()) {
            Atom atom = iterator.nextAtom();
            atoms[k++] = atom;
        } 
        for(int i=0; i<N-1; i++) {
            aPairs[i] = new AtomPair[N-1-i];
            for(int j=0; j<N-1-i; j++) {
                AtomPair aPair = new AtomPair();
                aPair.atom0 = atoms[i];
                aPair.atom1 = atoms[i+j+1];
                aPairs[i][j] = aPair;
            }
        }
    }
    
    private final AtomPair[][] aPairs;
}
