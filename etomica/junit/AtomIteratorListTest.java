package etomica.junit;

import java.util.Collections;
import java.util.LinkedList;

import junit.framework.TestCase;
import etomica.Atom;
import etomica.IteratorDirective;
import etomica.Phase;
import etomica.Simulation;
import etomica.Space;
import etomica.Species;
import etomica.SpeciesSpheresMono;
import etomica.atom.AtomLinker;
import etomica.atom.AtomList;
import etomica.atom.iterator.AtomIteratorList;
import etomica.space3d.Space3D;

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
	private Lister[] lister;
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
		Space space = new Space3D();
		Simulation sim = new Simulation();
		phase = new Phase();
		species = new SpeciesSpheresMono();
		species.setNMolecules(nMolecules);;
		sim.elementCoordinator.go();
		basis1 = species.getAgent(phase);
//		atomList = phase.speciesMaster.atomList;
		atomList = new AtomList();
		for(int i=0; i<nMolecules; i++) atomList.add(new Atom(space));
		first = atomList.getFirst();
		last = atomList.getLast();
		mid=atomList.size()/2;
		middle = atomList.get(mid);
		iterator = new AtomIteratorList(atomList);
		header = atomList.header;
		firstLinker = atomList.firstEntry();
		lastLinker = atomList.lastEntry();
		midLinker = atomList.entry(mid);
		lister=new Lister[nLists];
		for (int i=0;i<lister.length;i++) {
			lister[i] = new Lister();
		}
		list = new LinkedList[6];
	}
	
/**
 * Clears all of the lists.
 */
	private void clearLists() {
		for (int i=0;i<nLists;i++) {
			lister[i].list.clear();
		}
	}
	
/**
 * Prints all of the lists.
 */
	private void printLists() {
		for (int i=0;i<nLists;i++) {
			System.out.println(lister[i].list);
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
		assertTrue(lister[0].list.containsAll(lister[2].list));

		//test whether iterator does same thing twice
		assertEquals(lister[0].list,lister[1].list);
		assertEquals(lister[2].list,lister[3].list);
		
		Collections.reverse(lister[2].list);
		assertEquals(lister[0].list,lister[2].list);
		
		Collections.reverse(lister[1].list);
		assertEquals(lister[1].list,lister[3].list);
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
		assertEquals(lister[1].list,lister[0].list);
		

		iterator.setFirst(middle);
		iterator.allAtoms(lister[2]);
		lister[0].list.subList(0,mid).clear();
		assertEquals(lister[2].list,lister[0].list);

		iterator.setFirst(last);
		iterator.allAtoms(lister[3]);
		LinkedList list3 = new LinkedList();
		list3.add(last.toString());
		assertEquals(lister[3].list,list3);

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
		assertEquals(lister[2].list,list2);

		iterator.setFirst(middle);
		iterator.allAtoms(lister[3]);
		if(lister[0].list.size()%2==0) {
			lister[0].list.subList(0,mid-1).clear();
		}
		else {
			lister[0].list.subList(0,mid).clear();
		}
		assertEquals(lister[3].list,lister[0].list);

		iterator.setFirst(last);
		iterator.allAtoms(lister[4]);
		assertEquals(lister[4].list,lister[1].list);
		
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

		assertEquals(lister[0].list,lister[1].list);
		
		//Tests while loop iteration with atom set to first.
		iterator.setFirst(first);
		iterator.reset();
		while (iterator.hasNext()) {
			lister[2].actionPerformed(iterator.next());
		}
		assertEquals(lister[0].list,lister[2].list);
		
		//Tests while loop iteration with atom set to middle
		iterator.setFirst(middle);
		iterator.reset();
		while (iterator.hasNext()) {
			lister[3].actionPerformed(iterator.next());
		}
		assertEquals(lister[0].list.subList(mid,nMolecules),lister[3].list);

		//Tests while loop iteration with atom set to last
		iterator.setFirst(last);
		iterator.reset();
		while (iterator.hasNext()) {
			lister[4].actionPerformed(iterator.next());
		}
		assertEquals(lister[0].list.subList(nMolecules-1,nMolecules),lister[4].list);

		//SkipFirstAtom does not have any affect on the iterator
		
//		//Tests while loop iteration with setSkipFirst true
//		id.set();
//		iterator.reset(id);
//		iterator.setSkipFirstAtom(true);
//		while (iterator.hasNext()) {
//			lister[5].actionPerformed(iterator.next());
//		}
//		System.out.println(list[5]);
//		assertEquals(lister[0].list.subList(1,nMolecules),list[5]);

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
		assertEquals(lister[0].list,lister[1].list);
		
		//Test while loop iteration with atom set to last
		iterator.setFirst(last);
		iterator.reset();
		while (iterator.hasNext()) {
			lister[2].actionPerformed(iterator.next());
		}
		assertEquals(lister[0].list,lister[2].list);
		
		//Test while loop iteration with atom set to middle
		iterator.setFirst(middle);
		iterator.reset();
		while (iterator.hasNext()) {
			lister[3].actionPerformed(iterator.next());
		}
		if (nMolecules%2==0) {
			assertEquals(lister[0].list.subList(mid-1,nMolecules),lister[3].list);
		}
		else {
			assertEquals(lister[0].list.subList(mid,nMolecules),lister[3].list);
		}
		
		//Tests while loop iteration with atom set to first
		iterator.setFirst(first);
		iterator.reset();
		while (iterator.hasNext()) {
			lister[4].actionPerformed(iterator.next());
		}
		assertEquals(lister[0].list.subList(nMolecules-1,nMolecules),lister[4].list);
		
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
		list[3] = IteratorTest.generalIteratorTest(iterator);
		
		printLists();
		
		clearLists();		
	}

}
