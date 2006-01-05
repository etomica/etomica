package etomica.virial;

import etomica.atom.Atom;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomPositionFirstAtom;
import etomica.atom.iterator.AtomIteratorArrayListSimple;
import etomica.space.CoordinatePair;
import etomica.space.Space;

/**
 * @author David Kofke
 *
 * Class that holds a set of atom pairs.  Takes a list of atoms in its
 * constructor, and forms an instance of CoordinatePair for each pair formed in the
 * list.  Each CoordinatePair instance can be accessed via the getPair method.  The
 * CoordinatePairSet has an ID, which changes whenever the atom positions change 
 * (after which reset() should be called).  Cluster.value() depends on a one-to-one
 * correspondence between this ID and the positions of atoms in the CoordinatePairSet. 
 */

/* History
 * 08/21/03 (DAK) rearranged loops in resetPairs method
 */
 
public class CoordinatePairSet implements java.io.Serializable {

    /**
     * Constructor for CoordinatePairSet.
     * @param list The list of atoms for which the set of pairs is formed.
     */
    public CoordinatePairSet(AtomArrayList list, Space space) {
        cPairs = new CoordinatePair[list.size()][];
        setAtoms(list, space);
    }
    
    /**
     * Returns atom pair for ith and jth atoms in set.
     */
    public CoordinatePair getCPair(int i, int j) {
        if(i==j) throw new IllegalArgumentException("Error: asking for pair formed with both atoms the same");
        return i<j ? cPairs[i][j-i-1] : cPairs[j][i-j-1];
    }

    private void setAtoms(AtomArrayList list, Space space) {
        AtomIteratorArrayListSimple iterator = new AtomIteratorArrayListSimple(list);
        int N = list.size();
        Atom[] atoms = new Atom[N];
        iterator.reset();
        int k=0;
        while(iterator.hasNext()) {
            Atom atom = iterator.nextAtom();
            atoms[k++] = atom;
        }
        for(int i=0; i<N; i++) {
            cPairs[i] = new CoordinatePair[N-1-i];
            for(int j=0; j<N-1-i; j++) {
                CoordinatePair cPair;
                if (atoms[i].coord != null) {
                    cPair = new CoordinatePair(space);
                    cPair.reset(atoms[i].coord,atoms[i+j+1].coord);
                }
                else {
                    cPair = new CoordinatePairMolecular(space,new AtomPositionFirstAtom());
                    ((CoordinatePairMolecular)cPair).reset(atoms[i],atoms[i+j+1]);
                }
                cPairs[i][j] = cPair;
            }
        }
    }
    
    public void reset() {
        if (dirtyAtom == -1) {
            int N = cPairs.length;
            for(int i=0; i<N; i++) {
                CoordinatePair[] cPairs2 = cPairs[i];
                for(int j=0; j<N-i-1; j++) {
                    cPairs2[j].reset();
                }
            }
        }
        else {
            // reset only pairs containing dirty atom
            for (int i=0; i<dirtyAtom; i++) {
                cPairs[i][dirtyAtom-i-1].reset();
            }
            final CoordinatePair[] cPairs2 = cPairs[dirtyAtom];
            for (int i=0; i<cPairs2.length; i++) {
                cPairs2[i].reset();
            }
        }
        ID = staticID++;
    }
    
    public void E(CoordinatePairSet c) {
        int N = cPairs.length;
        for(int i=0; i<N; i++) {
            CoordinatePair[] cPairs2 = cPairs[i];
            CoordinatePair[] c2 = c.cPairs[i];
            for(int j=0; j<N-i-1; j++) {
                cPairs2[j].dr().E(c2[j].dr());
            }
        }
    }
    
    public int getID() {
        return ID;
    }
    
    protected final CoordinatePair[][] cPairs;
    private int ID;
    private static int staticID;
    public int dirtyAtom;
}