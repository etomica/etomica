package etomica.virial;

import etomica.*;

/**
 * @author David Kofke
 *
 * Class that holds a set of atom pairs.  Takes a list of atoms in its
 * constructor, and forms an instance of AtomPair for each pair formed in the
 * list.  Each AtomPair instance can be accessed via the getPair method, which
 * takes as arguments the node.index values for the two atoms of interest.
 */
public class PairSet {

	/**
	 * Constructor for PairSet.
	 * @param list The list of atoms for which the set of pairs is formed.
	 * Values of node.index for atoms in list must be in sequence with the
	 * atoms, i.e. the first atom has node.index = 0, the second has node.index
	 * = 1, etc.
	 */
	public PairSet(AtomList list) {
		super();
		setAtoms(list);
	}
	
	/**
	 * Returns atom pair for atoms having index i and j, respectively.
	 */
	public AtomPair getPair(int i, int j) {
		if(i==j) throw new IllegalArgumentException("Error: asking for pair formed with both atoms the same");
		return (i<j) ? pairs[i][j-i-1] : pairs[j][i-j-1];
	}

	private void setAtoms(AtomList list) {
		AtomIteratorList iterator = new AtomIteratorList(list);
		N = list.size();
		Atom[] atoms = new Atom[N];
		pairs = new AtomPair[N-1][];
		iterator.reset();
		int k=0;
		while(iterator.hasNext()) {
			Atom atom = iterator.next();
			if(atom.node.index() != k) throw new RuntimeException("Error in indexing of atoms");
			atoms[k++] = atom;
		} 
		for(int i=0; i<N-1; i++) {
			pairs[i] = new AtomPair[N-1-i];
			for(int j=0; j<N-1-i; j++) {
				pairs[i][j] = new AtomPair(atoms[i], atoms[i+j+1]);} 
		}
	}//end of setAtoms
	
	/**
	 * Calls  reset method for all pairs, computing all r-squared values for
	 * current configuration.  Returns this pairSet, so reset can be called in-
	 * line.
	 */
	public PairSet resetPairs() {
		for(int i=0; i<N-1; i++) {
			for(int j=0; j<N-1-i; j++) pairs[i][j].reset();
		}
		return this;
	}
	
	public int nMolecules() {return N;}

	private AtomPair[][] pairs;
	private int N;

}//end of PairSet
