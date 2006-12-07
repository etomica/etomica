package etomica.virial;

import etomica.atom.Atom;
import etomica.atom.AtomArrayList;
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
 
public class CoordinatePairMoleculeSet implements java.io.Serializable, CoordinatePairSet {

    /**
     * Constructor for CoordinatePairLeafSet.
     * @param list The list of atoms for which the set of pairs is formed.
     */
    public CoordinatePairMoleculeSet(AtomArrayList list, Space space) {
        atoms = new Atom[list.size()];
        numAtoms = list.size();
        r2 = new double[numAtoms*numAtoms];
        setAtoms(list);
        dr = space.makeVector();
        iPosition = space.makeVector();
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
            atoms[k++] = iterator.nextAtom();
        }
    }
    
    public void reset() {
        for(int i=0; i<numAtoms-1; i++) {
            Atom iAtom = atoms[i];
            iPosition.E(iAtom.getType().getPositionDefinition().position(iAtom));
            for(int j=i+1; j<numAtoms; j++) {
                Atom jAtom = atoms[j];
                Vector jPosition = jAtom.getType().getPositionDefinition().position(jAtom);
                dr.Ev1Mv2(iPosition, jPosition);
                r2[i*numAtoms+j] = dr.squared();
            }
        }
        ID = staticID++;
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
    protected final Atom[] atoms;
    protected final int numAtoms;
    protected final Vector dr;
    private final Vector iPosition;
    private int ID;
    private static int staticID;
}