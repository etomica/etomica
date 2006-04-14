package etomica.virial;

import etomica.atom.Atom;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomLeaf;
import etomica.atom.iterator.AtomIteratorArrayListSimple;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.util.Debug;

/**
 * @author David Kofke
 *
 * Class that holds a set of atom pairs.  Takes a list of atoms in its
 * constructor, and forms an instance of CoordinatePair for each pair formed in the
 * list.  Each CoordinatePair instance can be accessed via the getPair method.  The
 * CoordinatePairLeafSet has an ID, which changes whenever the atom positions change 
 * (after which reset() should be called).  Cluster.value() depends on a one-to-one
 * correspondence between this ID and the positions of atoms in the CoordinatePairLeafSet. 
 */

/* History
 * 08/21/03 (DAK) rearranged loops in resetPairs method
 */
 
public class CoordinatePairLeafSet implements java.io.Serializable, CoordinatePairSet {

    /**
     * Constructor for CoordinatePairLeafSet.
     * @param list The list of atoms for which the set of pairs is formed.
     */
    public CoordinatePairLeafSet(AtomArrayList list, Space space) {
        positions = new Vector[list.size()];
        numAtoms = list.size();
        r2 = new double[numAtoms*(numAtoms-1)/2];
        setAtoms(list);
        dirtyAtom = -1;
        dr = space.makeVector();
    }
    
    /**
     * Returns atom pair for ith and jth atoms in set.
     */
    public double getr2(int i, int j) {
        if(Debug.ON && !(i<j)) throw new IllegalArgumentException("Error: i must be less than j");
        return r2[i*numAtoms+j];
    }

    private void setAtoms(AtomArrayList list) {
        AtomIteratorArrayListSimple iterator = new AtomIteratorArrayListSimple(list);
        iterator.reset();
        int k=0;
        while(iterator.hasNext()) {
            Atom atom = iterator.nextAtom();
            positions[k++] = ((AtomLeaf)atom).coord.position();
        }
    }
    
    public void reset() {
        if (dirtyAtom == -1) {
            for(int i=0; i<numAtoms-1; i++) {
                Vector pos1 = positions[i];
                for(int j=i+1; j<numAtoms; j++) {
                    dr.Ev1Mv2(pos1,positions[j]);
                    r2[i*numAtoms+j] = dr.squared();
                }
            }
        }
        else {
            // reset only pairs containing dirty atom
            Vector posi = positions[dirtyAtom];
            for (int i=0; i<dirtyAtom; i++) {
                dr.Ev1Mv2(posi,positions[i]);
                r2[i*numAtoms+dirtyAtom] = dr.squared();
            }
            for (int j=dirtyAtom+1; j<numAtoms; j++) {
                dr.Ev1Mv2(posi,positions[j]);
                r2[dirtyAtom*numAtoms+j] = dr.squared();
            }
        }
        ID = staticID++;
    }
    
    public void setDirtyAtom(int newDirtyAtom) {
        dirtyAtom = newDirtyAtom;
    }
    
    public void E(CoordinatePairLeafSet c) {
        for(int i=0; i<numAtoms-1; i++) {
            for(int j=i+1; j<numAtoms; j++) {
                r2[i*numAtoms+j] =  c.r2[i*numAtoms+j];
            }
        }
    }
    
    public int getID() {
        return ID;
    }
    
    protected final double[] r2;
    protected final Vector[] positions;
    protected final int numAtoms;
    protected final Vector dr;
    private int ID;
    private static int staticID;
    protected int dirtyAtom;
}