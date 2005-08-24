package etomica.config;

import etomica.atom.Atom;
import etomica.atom.AtomList;
import etomica.atom.iterator.AtomIteratorListSimple;
import etomica.space.Space;
import etomica.space.Vector;

/**
 * General class for a collection of linearly linked atoms.
 * 
 * @author nancycribbin
 */

public abstract class ConformationChain extends Conformation {

	public ConformationChain(Space space){	
		super(space);		
		atomIterator = new AtomIteratorListSimple();
		//orientationVector = space.makeVector();
		//wrongNumberOfVectors = "Wrong number of vectors in the argument to ConformationChain subclass.";
	}
	
	/**
	 * Resets the vector instructions.
	 */
	protected abstract void reset();
	
	/**
	 * 
	 * @return the instructions to get to the location of the next molecule from
	 * the current one.
	 */
	protected abstract Vector nextVector();
	
	/**
	 * Places a set of atoms in a linearly connected fashion.
	 * @ param atomlist a list of atoms in the order in which they are linked
	 */
	public void initializePositions(AtomList atomlist){
		
		//First, check that we actually have some atoms
		int size = atomlist.size();
        	if(size == 0) return;
        
        	atomIterator.setList(atomlist);
        	atomIterator.reset();
        
        	reset();
		
        	//space.makeVector() zeroes the made Vector automatically
        	Vector currentPosition = space.makeVector();
        
        	//Zero the first atom.
        	atomIterator.nextAtom().coord.position().E(0.0);
        	
        	while(atomIterator.hasNext()){
        		Atom a = atomIterator.nextAtom();
        		//TODO someday, we might want a to be a chunk-of-atoms
        		currentPosition.PE(nextVector());
        		a.coord.position().E(currentPosition);
        	}
	}
	
	/**
	 * The vector drawn from the head of the molecule to the tail of the molecule.
	 */
/*	private Vector calculateOrientationVector(AtomList atomlist){
		//atomIterator.reset();
		Atom a = atomlist.getFirst();
		Atom z = atomlist.getLast();
		
		Vector ov = space.makeVector();
		
		ov.E(z.coord.position());
		ov.ME(a.coord.position());
		return ov;
	}
*/	
	
	//TODO Should we have a method here that moves the atom somehow- rotate, translate, etc.?
	
	private final AtomIteratorListSimple atomIterator;

	//private Vector orientationVector;	
}
