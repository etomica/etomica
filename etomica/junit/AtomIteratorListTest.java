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
	private AtomActive[] lister;
	private Atom basis1;
	private AtomList atomList;
	private AtomLinker.Tab header;
	private AtomLinker firstLinker,lastLinker,midLinker;
	private Atom first,middle,last;
	private int mid;

	public static void main(String[] args) {
		junit.textui.TestRunner.run(AtomIteratorListTest.class);
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
		lister=new AtomActive[nLists];
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
		iterator.all(basis1,new IteratorDirective(),lister[0]);
		iterator.all(basis1,new IteratorDirective(IteratorDirective.UP),lister[1]);
		iterator.all(basis1,new IteratorDirective(IteratorDirective.DOWN),lister[2]);
		iterator.all(basis1,new IteratorDirective(IteratorDirective.DOWN),lister[3]);
		iterator.all(basis1,new IteratorDirective(IteratorDirective.BOTH),lister[4]);
		iterator.all(basis1,new IteratorDirective(IteratorDirective.BOTH),lister[5]);
		iterator.all(basis1,new IteratorDirective(IteratorDirective.NEITHER),lister[6]);
		iterator.all(basis1,new IteratorDirective(IteratorDirective.NEITHER),lister[7]);

//		printLists();
		
		//test whether lists contain same elements
		assertTrue(list[0].containsAll(list[2]));
		assertTrue(list[2].containsAll(list[4]));

		//test whether iterator does same thing twice
		assertEquals(list[0],list[1]);
		assertEquals(list[2],list[3]);
		assertEquals(list[4],list[5]);
		assertEquals(list[6],list[7]);
		
		assertEquals(list[0],list[4]);
		
		Collections.reverse(list[2]);
		assertEquals(list[0],list[2]);
		
		Collections.reverse(list[1]);
		assertEquals(list[1],list[3]);
		assertEquals(list[6],new LinkedList());
		clearLists();
	}
	
/**
 * Tests the all method with UP Directives.
 */
	public void testUpDirectives() {

		IteratorDirective id = new IteratorDirective(IteratorDirective.UP);

		iterator.all(basis1,id,lister[0]);

		id.set(first);
		iterator.all(basis1,id,lister[1]);
		assertEquals(list[1],list[0]);

		id.set(middle);
		iterator.all(basis1,id,lister[2]);
		list[0].subList(0,mid).clear();
		assertEquals(list[2],list[0]);

		id.set(last);
		iterator.all(basis1,id,lister[3]);
		LinkedList list3 = new LinkedList();
		list3.add(last.toString());
		assertEquals(list[3],list3);
		
//		printLists();

		clearLists();
	}	

/**
 * Tests the all method with DOWN Directives
 */
	public void testDownDirectives() {

		IteratorDirective id = new IteratorDirective(IteratorDirective.DOWN);

		iterator.all(basis1,id,lister[0]);
		iterator.all(basis1,id,lister[1]);

		id.set(first);
		iterator.all(basis1,id,lister[2]);
		LinkedList list2 = new LinkedList();
		list2.add(first.toString());
		assertEquals(list[2],list2);

		id.set(middle);
		iterator.all(basis1,id,lister[3]);
		if(list[0].size()%2==0) {
		list[0].subList(0,mid-1).clear();
		}
		else {
		list[0].subList(0,mid).clear();
		}
		assertEquals(list[3],list[0]);

		id.set(last);
		iterator.all(basis1,id,lister[4]);
		assertEquals(list[4],list[1]);
		
//		printLists();
		
		clearLists();
	}
	
/**
 * Tests the all method with BOTH Directives.
 */
	public void testBothDirectives() {
		
		IteratorDirective id = new IteratorDirective(IteratorDirective.BOTH);

		iterator.all(basis1,id,lister[0]);
		iterator.all(basis1,id,lister[1]);

		id.set(first);
		iterator.all(basis1,id,lister[2]);
		assertEquals(list[2],list[0]);
		
		id.set(middle);
		iterator.all(basis1,id,lister[3]);
		Collections.reverse(list[0].subList(0,mid));
		Collections.rotate(list[0],-mid);
		assertEquals(list[0],list[3]);
		
		id.set(last);
		iterator.all(basis1,id,lister[4]);
		Collections.reverse(list[1].subList(0,list[1].size()-1));
		Collections.rotate(list[1],1);		
		assertEquals(list[4],list[1]);
		
//		printLists();
		
		clearLists();
	}
	
/**
 * Tests the all method with NEITHER Directives
 */
	public void testNeitherDirectives() {
		
		IteratorDirective id = new IteratorDirective(IteratorDirective.NEITHER);

		iterator.all(basis1,id,lister[0]);
		assertEquals(list[0],new LinkedList());
		
		id.set(first);
		iterator.all(basis1,id,lister[1]);
		LinkedList list1 = new LinkedList();
		list1.add(first.toString());
		assertEquals(list[1],list1);
		
		clearLists();
	}
	
/**
 * Tests a next loop with an up Directive.
 */
	public void testUpNextLoop() {

		IteratorDirective id = new IteratorDirective(IteratorDirective.UP);

		iterator.setBasis(basis1);	
		iterator.reset(id);
		iterator.all(basis1,id,lister[0]);
				
		//Tests Basic while loop iteration
		while (iterator.hasNext()) {
			lister[1].actionPerformed(iterator.next());
		}

		assertEquals(list[0],list[1]);
		
		//Tests while loop iteration with atom set to first.
		id.set(first);
		iterator.reset(id);
		while (iterator.hasNext()) {
			lister[2].actionPerformed(iterator.next());
		}
		assertEquals(list[0],list[2]);
		
		//Tests while loop iteration with atom set to middle
		id.set(middle);
		iterator.reset(id);
		while (iterator.hasNext()) {
			lister[3].actionPerformed(iterator.next());
		}
		assertEquals(list[0].subList(mid,nMolecules),list[3]);

		//Tests while loop iteration with atom set to last
		id.set(last);
		iterator.reset(id);
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
		
		IteratorDirective id = new IteratorDirective(IteratorDirective.DOWN);
		
		iterator.all(basis1,id,lister[0]);
		
		//Tests basic while loop iteration
		iterator.reset(id);
		while (iterator.hasNext()) {
			lister[1].actionPerformed(iterator.next());
		}
		assertEquals(list[0],list[1]);
		
		//Test while loop iteration with atom set to last
		id.set(last);
		iterator.reset(id);
		while (iterator.hasNext()) {
			lister[2].actionPerformed(iterator.next());
		}
		assertEquals(list[0],list[2]);
		
		//Test while loop iteration with atom set to middle
		id.set(middle);
		iterator.reset(id);
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
		id.set(first);
		iterator.reset(id);
		while (iterator.hasNext()) {
			lister[4].actionPerformed(iterator.next());
		}
		assertEquals(list[0].subList(nMolecules-1,nMolecules),list[4]);
		
		clearLists();
	}

/**
 * Tests a next loop with a BOTH Directive
 */
	public void testBothNextLoop() {
		
		IteratorDirective id = new IteratorDirective(IteratorDirective.BOTH);
		
		iterator.reset(id);
		
		iterator.all(basis1,id,lister[0]);
		
		//Tests basic while loop iteration
		while (iterator.hasNext()) {
			lister[1].actionPerformed(iterator.next());
		}
		assertEquals(list[0],list[1]);
		
		//Tests while loop iteration with atom set to first
		id.set(first);
		iterator.reset(id);
		while (iterator.hasNext()) {
			lister[2].actionPerformed(iterator.next());
		}
		assertEquals(list[0],list[1]);
		
		//Tests while loop iteration with atom set to last
		id.set(last);
		iterator.reset(id);
		while (iterator.hasNext()) {
			lister[3].actionPerformed(iterator.next());
		}
		iterator.all(basis1,new IteratorDirective(IteratorDirective.DOWN),lister[4]);
		assertEquals(list[4],list[3]);
		
		//Tests while loop iteration with atom set to middle
		id.set(middle);
		iterator.reset(id);
		while (iterator.hasNext()) {
			lister[5].actionPerformed(iterator.next());
		}
		Collections.reverse(list[0].subList(0,mid));
		Collections.rotate(list[0],-mid);
		assertEquals(list[0],list[5]);
		
		clearLists();
	}
	
/**
 * Tests a next loop with a NEITHER Directive
 */
	public void testNeitherNextLoop() {
		
		IteratorDirective id = new IteratorDirective(IteratorDirective.NEITHER);
		
		iterator.reset(id);
		
		assertEquals(false,iterator.hasNext());
		
		id.set(middle);
		iterator.reset(id);
		
		assertEquals(middle,iterator.next());
		assertEquals(false,iterator.hasNext());
	}
	
/**
 * Tests the other all Methods with an UP Directive
 */
	public void testOtherAllMethodsUP() {
		
		IteratorDirective id = new IteratorDirective();
		
		//Tests from first atom/linker
		iterator.all(basis1,id,lister[0]);
		AtomIteratorList.all(firstLinker,null,id,lister[1]);
		AtomIteratorList.all(firstLinker,id,lister[2]);
		AtomIteratorList.all(firstLinker,header,id,lister[3]);
//		printLists();
		assertEquals(list[0],list[1]);
		assertEquals(list[1],list[2]);
		assertEquals(list[2],list[3]);
		
		clearLists();
		
		//Tests from middle atom/linker
		id.set(middle);
		
		iterator.all(basis1,id,lister[0]);
		AtomIteratorList.all(midLinker,null,id,lister[1]);
		AtomIteratorList.all(midLinker,id,lister[2]);
		AtomIteratorList.all(midLinker,header,id,lister[3]);
		assertEquals(list[0],list[1]);
		assertEquals(list[1],list[2]);
		assertEquals(list[2],list[3]);
				
		clearLists();
		
		//Test from last atom/linker
		id.set(last);
		
		iterator.all(basis1,id,lister[0]);
		AtomIteratorList.all(lastLinker,null,id,lister[1]);
		AtomIteratorList.all(lastLinker,id,lister[2]);
		AtomIteratorList.all(lastLinker,header,id,lister[3]);
		assertEquals(list[0],list[1]);
		assertEquals(list[1],list[2]);
		assertEquals(list[2],list[3]);		
		
		clearLists();
	}
	
/**
 * Tests the other all methods with a DOWN directive
 */
	public void testOtherAllMethodsDOWN() {
		
		IteratorDirective id = new IteratorDirective(IteratorDirective.DOWN);
		
		//Tests from last atom/linker
		iterator.all(basis1,id,lister[0]);
		AtomIteratorList.all(lastLinker,null,id,lister[1]);
		AtomIteratorList.all(lastLinker,id,lister[2]);
		AtomIteratorList.all(lastLinker,header,id,lister[3]);
		assertEquals(list[0],list[1]);
		assertEquals(list[1],list[2]);
		assertEquals(list[2],list[3]);
		
		clearLists();
		
		//Tests from middle atom/linker
		id.set(middle);
		
		iterator.all(basis1,id,lister[0]);
		AtomIteratorList.all(midLinker,null,id,lister[1]);
		AtomIteratorList.all(midLinker,id,lister[2]);
		AtomIteratorList.all(midLinker,header,id,lister[3]);
		assertEquals(list[0],list[1]);
		assertEquals(list[1],list[2]);
		assertEquals(list[2],list[3]);
		
		clearLists();
		
		//Test from first atom/linker
		id.set(first);
		
		iterator.all(basis1,id,lister[0]);
		AtomIteratorList.all(firstLinker,null,id,lister[1]);
		AtomIteratorList.all(firstLinker,id,lister[2]);
		AtomIteratorList.all(firstLinker,header,id,lister[3]);
		assertEquals(list[0],list[1]);
		assertEquals(list[1],list[2]);
		assertEquals(list[2],list[3]);
		
		clearLists();
	}
	
/**
 * Tests other all methods with a BOTH directive
 */
	public void testOtherAllMethodsBOTH() {
		
		IteratorDirective id = new IteratorDirective(IteratorDirective.BOTH);

		//Tests from first atom/linker		
		iterator.all(basis1,id,lister[0]);
		AtomIteratorList.all(firstLinker,null,id,lister[1]);
		AtomIteratorList.all(firstLinker,id,lister[2]);
		AtomIteratorList.all(firstLinker,header,id,lister[3]);
		assertEquals(list[0],list[1]);
		assertEquals(list[1],list[2]);
		assertEquals(list[2],list[3]);
		
		clearLists();

		//Tests from middle atom/linker
		id.set(middle);
		
		iterator.all(basis1,id,lister[0]);
		AtomIteratorList.all(midLinker,null,id,lister[1]);
		AtomIteratorList.all(midLinker,id,lister[2]);
		AtomIteratorList.all(midLinker,header,id,lister[3]);
		assertEquals(list[0],list[1]);
		assertEquals(list[1],list[2]);
		assertEquals(list[2],list[3]);
		
		clearLists();
		
		//Tests from last atom/linker
		id.set(last);
		
		iterator.all(basis1,id,lister[0]);
		AtomIteratorList.all(lastLinker,null,id,lister[1]);
		AtomIteratorList.all(lastLinker,id,lister[2]);
		AtomIteratorList.all(lastLinker,header,id,lister[3]);
		assertEquals(list[0],list[1]);
		assertEquals(list[1],list[2]);
		assertEquals(list[2],list[3]);
		
		clearLists();
	}
	
/**
 * Tests other all methods with a NEITHER Directive
 */
	public void testOtherAllMethdosNEITHER() {
		
		IteratorDirective id = new IteratorDirective(IteratorDirective.NEITHER);
		id.set(first);

		iterator.all(basis1,id,lister[0]);
		AtomIteratorList.all(firstLinker,null,id,lister[1]);
		AtomIteratorList.all(firstLinker,id,lister[2]);
		AtomIteratorList.all(firstLinker,header,id,lister[3]);
		
		assertEquals(list[0],list[1]);
		assertEquals(list[1],list[2]);
		assertEquals(list[2],list[3]);
		
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
		iterator.reset(midLinker,header,IteratorDirective.NEITHER);
		iterator.reset();
		while (iterator.hasNext()) {
			lister[1].actionPerformed(iterator.next());
		}
		assertEquals(list[0],list[1]);
		
		//tests the reset(int) method
		iterator.reset();
		iterator.reset(mid);
		while (iterator.hasNext()) {
			lister[2].actionPerformed(iterator.next());
		}
		AtomIteratorList.all(midLinker,new IteratorDirective(),lister[3]);
		assertEquals(list[2],list[3]);
		
		//tests the reset(Atom) method
		iterator.reset();
		iterator.reset(middle);
		while (iterator.hasNext()) {
			lister[4].actionPerformed(iterator.next());
		}
		AtomIteratorList.all(midLinker,new IteratorDirective(),lister[5]);
		assertEquals(list[4],list[5]);

		//tests the reset(AtomLinker) method
		iterator.reset();
		iterator.reset(midLinker);
		while (iterator.hasNext()) {
			lister[6].actionPerformed(iterator.next());
		}
		AtomIteratorList.all(midLinker,new IteratorDirective(),lister[7]);
		assertEquals(list[6],list[7]);
		
		clearLists();
				
		//tests the reset(Atom,Direction) method
		iterator.reset();
		iterator.reset(middle,IteratorDirective.DOWN);
		while (iterator.hasNext()) {
			lister[0].actionPerformed(iterator.next());
		}
		AtomIteratorList.all(midLinker,new IteratorDirective(IteratorDirective.DOWN),lister[1]);
		assertEquals(list[0],list[1]);
		
		//tests the reset(AtomLinker,Direction) method
		iterator.reset();
		iterator.reset(midLinker,IteratorDirective.DOWN);
		while (iterator.hasNext()) {
			lister[2].actionPerformed(iterator.next());
		}
		AtomIteratorList.all(midLinker,new IteratorDirective(IteratorDirective.DOWN),lister[3]);
		assertEquals(list[2],list[3]);
		
		//tests the reset(AtomLinker,AtomLinker.Tab) method
		iterator.reset();
		iterator.reset(midLinker,header);
		while (iterator.hasNext()) {
			lister[4].actionPerformed(iterator.next());
		}
		AtomIteratorList.all(midLinker,new IteratorDirective(),lister[5]);
		assertEquals(list[4],list[5]);
		
		//tests the reset(AtomLinker,AtomLiner.Tab,Direction) method
		iterator.reset();
		iterator.reset(midLinker,header,IteratorDirective.DOWN);
		while (iterator.hasNext()) {
			lister[6].actionPerformed(iterator.next());
		}
		AtomIteratorList.all(midLinker,new IteratorDirective(IteratorDirective.DOWN),lister[7]);
		assertEquals(list[6],list[7]);
		
		clearLists();		
	}
	
	public void testUnset() {
		iterator.unset();
		assertEquals(false,iterator.hasNext());
	}

/**
 * 
 * @author Andrew Walker
 *
 * This class creates an AtomActive which puts the result of an Atom's toString
 * method into a list in the array of lists in the Testing Class
 */
	class Lister implements AtomActive {

		private final int i;

		public Lister(int i) {
			this.i=i;
		}
		public void actionPerformed(AtomSet as) {
			list[i].add(as.toString());	
		}
		public void actionPerformed(Atom a) {
			list[i].add(a.toString());
		}
	}

}
