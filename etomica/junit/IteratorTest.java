package etomica.junit;

import java.util.Collections;
import java.util.LinkedList;

import etomica.*;
import junit.framework.*;

/**
 * @author aawalker
 *
 */
abstract class IteratorTest extends TestCase {
	
	protected Simulation sim;
	protected Lister[] lister;
	protected Species species1,species2,species3;
	protected int nMolecules=15;
	protected int nAtoms=3;
	protected int nLists=8; //minimum 8
	protected Phase phase;
	protected Atom basis1,basis2;
	protected AtomPair pairBasis;
	protected AtomList atomList,atomList1,atomList2;
	protected AtomLinker.Tab header;
	protected AtomLinker firstLinker1,lastLinker1,midLinker1,firstLinker2,lastLinker2,midLinker2;
	protected Atom first1,middle1,last1,first2,middle2,last2,outsideAtom;
	protected int mid1,mid2;
	
	/**
	 * This method is called before any test is run and sets up the basics of the
	 * etomica simulation and testing elements.
	 */
	protected void setUp() {
		sim = new Simulation();
		phase = new Phase();
		species1 = new SpeciesSpheresMono(nMolecules);
		species2 = new SpeciesSpheresMono(nMolecules);
		species3 = new SpeciesSpheresMono(1);
		sim.elementCoordinator.go();
		lister = Lister.listerArray(nLists);
		basis1 = species1.getAgent(phase);
		basis2 = species2.getAgent(phase);
		pairBasis = new AtomPair(basis1,basis2);
		atomList = phase.speciesMaster.atomList;
		atomList1 = ((AtomTreeNodeGroup)basis1.node).childList;
		atomList2 = ((AtomTreeNodeGroup)basis2.node).childList;
		first1 = atomList1.getFirst();
		last1 = atomList1.getLast();
		mid1=atomList1.size()/2;
		middle1 = atomList1.get(mid1);
		first2 = atomList2.getFirst();
		last2 = atomList2.getLast();
		mid2 = atomList2.size()/2;
		middle2 = atomList2.get(mid2);
		header = atomList.header;
		firstLinker1 = atomList1.firstEntry();
		lastLinker1 = atomList1.lastEntry();
		midLinker1 = atomList1.entry(mid1);
		firstLinker2 = atomList2.firstEntry();
		lastLinker2 = atomList2.lastEntry();
		midLinker2 = atomList2.entry(mid2);
		outsideAtom = species3.getAgent(phase).firstMolecule();
	}
	
	/**
	 * Clears all of the lists.
	 */
		protected void clearLists() {
			for (int i=0;i<nLists;i++) {
				lister[i].list.clear();
			}
		}
	
	/**
	 * Prints all of the lists.
	 */
		protected void printLists() {
			for (int i=0;i<nLists;i++) {
				System.out.println(lister[i].list);
			}
			System.out.println();
		}
		
		public static java.util.LinkedList generalIteratorTest(AtomsetIterator iterator) {
			Lister[] lister = Lister.listerArray(4);
			iterator.allAtoms(lister[0]);
			
			AtomIterator atomIterator = (iterator instanceof AtomIterator) ? (AtomIterator)iterator : null;
			
			iterator.reset();
			while(iterator.hasNext()) {
				Atom[] peekAtom = iterator.peek();
				assertTrue(java.util.Arrays.equals(peekAtom, iterator.next()));
				lister[1].actionPerformed(peekAtom);
			}
			
			//test that allAtoms and hasNext/next give same set of iterates
			assertEquals(lister[0].list, lister[1].list);
			
			//test operation of unset method
			iterator.reset();
			iterator.unset();
			assertFalse(iterator.hasNext());
			assertNull(iterator.next()[0]);
			assertFalse(iterator.hasNext());
			if(atomIterator != null) assertNull(atomIterator.nextAtom());
			assertFalse(iterator.hasNext());
			
			//test size method
			assertEquals(iterator.size(), lister[0].list.size());
			
			//test contains method
			if(atomIterator != null) {
				AtomList atomList = new AtomList(atomIterator);
				AtomIteratorList listIterator = new AtomIteratorList(atomList);
				listIterator.reset();
				while(listIterator.hasNext()) {
					assertTrue(iterator.contains(listIterator.next()));
				}
				assertFalse(iterator.contains(null));
			}
			
			//test nBody
			iterator.reset();
			assertEquals(iterator.next().length, iterator.nBody());
						
			return lister[0].list;
		}
}
