package etomica.junit;

import etomica.*;
import junit.framework.*;
import java.util.*;

/**
 * @author Andrew Walker
 *
 * This class tests methods of AtomIteratorList using JUnit tests.  It should be
 * run from Eclipse by choosing "Run As > JUnit Test".
 */
public class AtomIteratorListTest extends TestCase {

	private Species species;
	private AtomIteratorList iterator;
	private final int nMolecules=15;
	private final int nLists=8; //minimum 8
	private Phase phase;
	private LinkedList[] list;
	private AtomsetActive[] lister;
	private Atom basis1;
	private AtomList atomList;
	private AtomLinker.Tab header;
	private AtomLinker firstLinker,lastLinker,midLinker;
	private Atom first,middle,last;
	private int mid;

	public static void main(String[] args) {
		AtomIteratorListTest test = new AtomIteratorListTest();
		test.setUp();
		Atom[] atoms = new Atom[3];
		System.out.println(atoms.toString());
//		junit.textui.TestRunner.run(AtomIteratorListTest.class);
	}

/**
 * This method is called before any test is run and sets up the basics of the
 * etomica simulation and testing elements.
 */
	protected void setUp() {
		Simulation sim = new Simulation();
		phase = new Phase();
		species = new SpeciesSpheresMono();
		species.setNMolecules(nMolecules);;
		sim.elementCoordinator.go();
		basis1 = species.getAgent(phase);
		atomList = phase.speciesMaster.atomList;
		first = atomList.getFirst();
		last = atomList.getLast();
		mid=atomList.size()/2;
		middle = atomList.get(mid);
		iterator = new AtomIteratorList(atomList);
		header = atomList.header;
		firstLinker = atomList.firstEntry();
		lastLinker = atomList.lastEntry();
		midLinker = atomList.entry(mid);
		list = new LinkedList[nLists];
		for (int i=0;i<nLists;i++) {
			list[i] = new LinkedList();
		}
		lister=new AtomsetActive[nLists];
		for (int i=0;i<list.length;i++) {
			lister[i] = new Lister(i);
		}
	}
	
/**
 * Clears all of the lists.
 */
	private void clearLists() {
		for (int i=0;i<nLists;i++) {
			list[i].clear();
		}
	}
	
/**
 * Prints all of the lists.
 */
	private void printLists() {
		for (int i=0;i<nLists;i++) {
			System.out.println(list[i]);
		}
		System.out.println();
	}

/**
 * Tests the size() method of the iterator.
 */
	public void testSize() {
		assertTrue(iterator.size()==nMolecules);
	}

/**
 * Tests the all() method with IteratorDirectives of all different directions.
 */
	public void testAllDirections() {
		iterator.setFirst((AtomLinker)null);

		iterator.setDirection(IteratorDirective.UP);
		iterator.allAtoms(lister[0]);
		iterator.setDirection(IteratorDirective.UP);
		iterator.allAtoms(lister[1]);
		iterator.setDirection(IteratorDirective.DOWN);
		iterator.allAtoms(lister[2]);
		iterator.setDirection(IteratorDirective.DOWN);
		iterator.allAtoms(lister[3]);

//		printLists();
		
		//test whether lists contain same elements
		assertTrue(list[0].containsAll(list[2]));

		//test whether iterator does same thing twice
		assertEquals(list[0],list[1]);
		assertEquals(list[2],list[3]);
		
		Collections.reverse(list[2]);
		assertEquals(list[0],list[2]);
		
		Collections.reverse(list[1]);
		assertEquals(list[1],list[3]);
		clearLists();
	}
	
/**
 * Tests the all method with UP Directives.
 */
	public void testUpDirectives() {

		iterator.setDirection(IteratorDirective.UP);
		iterator.setFirst((AtomLinker)null);

		iterator.allAtoms(lister[0]);

		iterator.setFirst(first);
		iterator.allAtoms(lister[1]);
//		printLists();
		assertEquals(list[1],list[0]);
		

		iterator.setFirst(middle);
		iterator.allAtoms(lister[2]);
		list[0].subList(0,mid).clear();
		assertEquals(list[2],list[0]);

		iterator.setFirst(last);
		iterator.allAtoms(lister[3]);
		LinkedList list3 = new LinkedList();
		list3.add(last.toString());
		assertEquals(list[3],list3);

		clearLists();
	}	

/**
 * Tests the all method with DOWN Directives
 */
	public void testDownDirectives() {

		iterator.setDirection(IteratorDirective.DOWN);
		iterator.setFirst((AtomLinker)null);

		iterator.allAtoms(lister[0]);
		iterator.allAtoms(lister[1]);

		iterator.setFirst(first);
		iterator.allAtoms(lister[2]);
		LinkedList list2 = new LinkedList();
		list2.add(first.toString());
		assertEquals(list[2],list2);

		iterator.setFirst(middle);
		iterator.allAtoms(lister[3]);
		if(list[0].size()%2==0) {
			list[0].subList(0,mid-1).clear();
		}
		else {
			list[0].subList(0,mid).clear();
		}
		assertEquals(list[3],list[0]);

		iterator.setFirst(last);
		iterator.allAtoms(lister[4]);
		assertEquals(list[4],list[1]);
		
//		printLists();
		
		clearLists();
	}
	
		
/**
 * Tests a next loop with an up Directive.
 */
	public void testUpNextLoop() {

		iterator.setDirection(IteratorDirective.UP);
		iterator.setFirst((AtomLinker)null);

		iterator.allAtoms(lister[0]);
				
		//Tests Basic while loop iteration
		iterator.reset();
		while (iterator.hasNext()) {
			lister[1].actionPerformed(iterator.next());
		}

		assertEquals(list[0],list[1]);
		
		//Tests while loop iteration with atom set to first.
		iterator.setFirst(first);
		iterator.reset();
		while (iterator.hasNext()) {
			lister[2].actionPerformed(iterator.next());
		}
		assertEquals(list[0],list[2]);
		
		//Tests while loop iteration with atom set to middle
		iterator.setFirst(middle);
		iterator.reset();
		while (iterator.hasNext()) {
			lister[3].actionPerformed(iterator.next());
		}
		assertEquals(list[0].subList(mid,nMolecules),list[3]);

		//Tests while loop iteration with atom set to last
		iterator.setFirst(last);
		iterator.reset();
		while (iterator.hasNext()) {
			lister[4].actionPerformed(iterator.next());
		}
		assertEquals(list[0].subList(nMolecules-1,nMolecules),list[4]);

		//SkipFirstAtom does not have any affect on the iterator
		
//		//Tests while loop iteration with setSkipFirst true
//		id.set();
//		iterator.reset(id);
//		iterator.setSkipFirstAtom(true);
//		while (iterator.hasNext()) {
//			lister[5].actionPerformed(iterator.next());
//		}
//		System.out.println(list[5]);
//		assertEquals(list[0].subList(1,nMolecules),list[5]);

		clearLists();
	}
	
/**
 * Tests a next loop with a Down directive
 */
	public void testDownNextLoop() {
		
		iterator.setDirection(IteratorDirective.DOWN);
		iterator.setFirst((AtomLinker)null);
		
		iterator.allAtoms(lister[0]);
		
		//Tests basic while loop iteration
		iterator.reset();
		while (iterator.hasNext()) {
			lister[1].actionPerformed(iterator.next());
		}
		assertEquals(list[0],list[1]);
		
		//Test while loop iteration with atom set to last
		iterator.setFirst(last);
		iterator.reset();
		while (iterator.hasNext()) {
			lister[2].actionPerformed(iterator.next());
		}
		assertEquals(list[0],list[2]);
		
		//Test while loop iteration with atom set to middle
		iterator.setFirst(middle);
		iterator.reset();
		while (iterator.hasNext()) {
			lister[3].actionPerformed(iterator.next());
		}
		if (nMolecules%2==0) {
			assertEquals(list[0].subList(mid-1,nMolecules),list[3]);
		}
		else {
			assertEquals(list[0].subList(mid,nMolecules),list[3]);
		}
		
		//Tests while loop iteration with atom set to first
		iterator.setFirst(first);
		iterator.reset();
		while (iterator.hasNext()) {
			lister[4].actionPerformed(iterator.next());
		}
		assertEquals(list[0].subList(nMolecules-1,nMolecules),list[4]);
		
		clearLists();
	}

		
	
	
	
/**
 * Tests the different resets
 */	
	public void testResets() {
		
		//tests the basic reset() for a has/next loop
		iterator.reset();
		while (iterator.hasNext()) {
			lister[0].actionPerformed(iterator.next());
		}
		
		//tests the reset(Atom) method
		iterator.setDirection(IteratorDirective.UP);
		iterator.setFirst(middle);
		list[2] = IteratorTest.generalIteratorTest(iterator);
		
		
		//tests the reset(AtomLinker) method
		iterator.setFirst(midLinker);
		list[6] = IteratorTest.generalIteratorTest(iterator);
		
		printLists();
		
		clearLists();		
	}
	
	public void testUnset() {
		iterator.unset();
		assertEquals(false,iterator.hasNext());
		iterator.next();
		assertEquals(false,iterator.hasNext());		
	}

/**
 * 
 * @author Andrew Walker
 *
 * This class creates an AtomActive which puts the result of an Atom's toString
 * method into a list in the array of lists in the Testing Class
 */
	class Lister implements AtomsetActive {

		private final int i;

		public Lister(int i) {
			this.i=i;
		}
		public void actionPerformed(Atom[] a) {
			list[i].add(a[0].toString());
		}
	}

}
