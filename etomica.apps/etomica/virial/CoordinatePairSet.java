package etomica.virial;

import etomica.Atom;
import etomica.Space;
import etomica.atom.AtomList;
import etomica.atom.iterator.AtomIteratorList;
import etomica.space.CoordinatePair;

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
 
public class CoordinatePairSet {

	/**
	 * Constructor for CoordinatePairSet.
     * @param list The list of atoms for which the set of pairs is formed.
	 */
	public CoordinatePairSet(AtomList list, Space space) {
		cPairs = new CoordinatePair[list.size()-1][];
		setAtoms(list, space);
	}
	
    /**
     * Returns atom pair for ith and jth atoms in set.
     */
    public CoordinatePair getCPair(int i, int j) {
        if(i==j) throw new IllegalArgumentException("Error: asking for pair formed with both atoms the same");
        return i<j ? cPairs[i][j-i-1] : cPairs[j][i-j-1];
    }

	private void setAtoms(AtomList list, Space space) {
		AtomIteratorList iterator = new AtomIteratorList(list);
		int N = list.size();
		Atom[] atoms = new Atom[N];
		iterator.reset();
		int k=0;
		while(iterator.hasNext()) {
			Atom atom = iterator.nextAtom();
			atoms[k++] = atom;
		} 
		for(int i=0; i<N-1; i++) {
			cPairs[i] = new CoordinatePair[N-1-i];
			for(int j=0; j<N-1-i; j++) {
				CoordinatePair cPair = space.makeCoordinatePair();
				cPair.reset(atoms[i].coord,atoms[i+j+1].coord);
				cPairs[i][j] = cPair;
			}
		}
	}
		
	public void reset() {
		int N = cPairs.length;
		for(int i=0; i<N; i++) {
			CoordinatePair[] cPairs2 = cPairs[i];
			for(int j=0; j<N-i; j++) {
				cPairs2[j].reset();
			}
		}
		ID = staticID++;
	}
	
	public int getID() {
		return ID;
	}
	
    private final CoordinatePair[][] cPairs;
	private int ID;
	private static int staticID;
}