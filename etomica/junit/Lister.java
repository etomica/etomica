package etomica.junit;

import etomica.*;
import java.util.LinkedList;

/**
 * @author aawalker
 *
 */
class Lister implements AtomActive, AtomPairActive {
	
	public final LinkedList list;
	private Space space;
	
	public Lister(Space spc) {
		list = new LinkedList();
		space = spc;
	}

	/**
	 * @see etomica.AtomPairActive#actionPerformed(etomica.AtomPair)
	 */
	public void actionPerformed(AtomPair pair) {
		list.add(pair.toString());
	}

	/**
	 * @see etomica.AtomPairActive#innerWrapper()
	 */
	public InnerWrapper innerWrapper() {
		return new AtomPairActive.InnerWrapper(this,new AtomPair(space));
	}

	/**
	 * @see etomica.AtomPairActive#outerWrapper()
	 */
	public OuterWrapper outerWrapper() {
		return new AtomPairActive.OuterWrapper(this);
	}

	/**
	 * @see etomica.AtomSetActive#actionPerformed(etomica.AtomSet)
	 */
	public void actionPerformed(AtomSet atomSet) {
		list.add(atomSet.toString());
	}
	
	/**
	 * Gives an array of listers each with their own list.
	 * 
	 * @param n the number of listers in the array
	 * @return AtomLister[]
	 */	
	public static Lister[] listerArray(int n,Space space) {
		Lister[] lister = new Lister[n];
		for (int i=0;i<n;i++) {
			lister[i] = new Lister(space);
		}
		return lister;
	}
	
	/**
	 * @see etomica.AtomActive#actionPerformed(etomica.Atom)
	 */
	public void actionPerformed(Atom atom) {
		list.add(atom.toString());
	}

}
