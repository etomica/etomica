/*
 * History
 * Created on Oct 13, 2004 by kofke
 */
package etomica.nbr;

import etomica.Atom;
import etomica.AtomIterator;
import etomica.AtomIteratorAtomDependent;
import etomica.AtomLinker;
import etomica.AtomsetActive;
import etomica.lattice.AbstractCell;

/**
 * @author kofke
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
public class AtomIteratorCell implements AtomIteratorAtomDependent {

	/**
	 * 
	 */
	public AtomIteratorCell() {
		super();
		// TODO Auto-generated constructor stub
	}

	/* (non-Javadoc)
	 * @see etomica.AtomIteratorAtomDependent#setAtom(etomica.Atom)
	 */
	public void setAtom(Atom atom) {
		// TODO find atom's cell

	}

	/* (non-Javadoc)
	 * @see etomica.AtomIterator#nextAtom()
	 */
	public Atom nextAtom() {
		// TODO Auto-generated method stub
		return null;
	}

	/* (non-Javadoc)
	 * @see etomica.AtomsetIterator#contains(etomica.Atom[])
	 */
	public boolean contains(Atom[] atom) {
		// TODO Auto-generated method stub
		return false;
	}

	/* (non-Javadoc)
	 * @see etomica.AtomsetIterator#hasNext()
	 */
	public boolean hasNext() {
		// TODO Auto-generated method stub
		return false;
	}

	/* (non-Javadoc)
	 * @see etomica.AtomsetIterator#reset()
	 */
	public void reset() {
		// TODO Auto-generated method stub

	}

	/* (non-Javadoc)
	 * @see etomica.AtomsetIterator#unset()
	 */
	public void unset() {
		// TODO Auto-generated method stub

	}

	/* (non-Javadoc)
	 * @see etomica.AtomsetIterator#next()
	 */
	public Atom[] next() {
		// TODO Auto-generated method stub
		return null;
	}

	/* (non-Javadoc)
	 * @see etomica.AtomsetIterator#peek()
	 */
	public Atom[] peek() {
		// TODO Auto-generated method stub
		return null;
	}

	public void setCell(AbstractCell referenceCell) {
		this.referenceCell = referenceCell;
	}
	
	/* (non-Javadoc)
	 * @see etomica.AtomsetIterator#allAtoms(etomica.AtomsetActive)
	 */
	public void allAtoms(AtomsetActive action) {
		AtomLinker.Tab header = referenceCell.neighborManager().neighbors().header;
		for(AtomLinker e=header.next; e!=header; e=e.next) {//loop over cells
			AtomLinker.Tab[] tabs = (AtomLinker.Tab[])e.atom.agents[0];
			AtomLinker next = tabs[tabIndex].next;
			while(next.atom != null) {//loop over atoms in cell
				action.actionPerformed(next.atom);
				next = next.next;
			}//end while
		}
	}

	/* (non-Javadoc)
	 * @see etomica.AtomsetIterator#size()
	 */
	public int size() {
		// TODO Auto-generated method stub
		return 0;
	}

	/* (non-Javadoc)
	 * @see etomica.AtomsetIterator#nBody()
	 */
	public int nBody() {
		return 1;
	}
	
	public void assignCells() {
		// TODO loop over atomIterator iterates and assign each atom to a cell
	}

	AtomIterator atomIterator;
	private AbstractCell referenceCell;
	
}
