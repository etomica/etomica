package etomica.junit;

import etomica.*;
import junit.framework.*;
import java.util.*;

/**
 * @author aawalker
 *
 * To change this generated comment edit the template variable "typecomment":
 * Window>Preferences>Java>Templates.
 * To enable and disable the creation of type comments go to
 * Window>Preferences>Java>Code Generation.
 */
public class AtomIteratorListTest extends TestCase {

	private Species species;
	private AtomIterator iterator;
	private final int nMolecules=15;
	private final int nLists=8; //minimum 8
	private Phase phase;
	private LinkedList[] list;
	private AtomActive[] lister;
	private Atom basis1;
	private AtomList atomList;
	private Atom first,middle,last;
	private int mid;

	public static void main(String[] args) {
		junit.textui.TestRunner.run(AtomIteratorListTest.class);
	}

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
		iterator=new AtomIteratorList(atomList);
		list=new LinkedList[nLists];
		for (int i=0;i<nLists;i++) {
			list[i]=new LinkedList();
		}
		lister=new AtomActive[nLists];
		for (int i=0;i<list.length;i++) {
			lister[i]=new Lister(i);
		}
	}
	
	private void clearLists() {
		for (int i=0;i<nLists;i++) {
			list[i].clear();
		}
	}
	
	private void printLists() {
		for (int i=0;i<nLists;i++) {
			System.out.println(list[i]);
		}
		System.out.println();
	}

	public void testNMolecules() {
		assertTrue(iterator.size()==nMolecules);
	}

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
		//
		assertEquals(list[0],list[4]);
		//
		Collections.reverse(list[2]);
		assertEquals(list[0],list[2]);
		Collections.reverse(list[1]);
		assertEquals(list[1],list[3]);
		//
		assertEquals(list[6],new LinkedList());
		clearLists();
	}
	
	public void testSkipFirst() {
		IteratorDirective id = new IteratorDirective();
		iterator.all(basis1,id,lister[0]);
		list[0].remove(0);
		id.setSkipFirst(true);
		iterator.all(basis1,id,lister[1]);
		assertEquals(list[0],list[1]);
		clearLists();
	}
	
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
	
	public void testClear() {
		IteratorDirective id = new IteratorDirective();
		id.set(last);
		id.set(IteratorDirective.DOWN);
		iterator.all(basis1,id,lister[0]);
		id.clear();
		iterator.all(basis1,id,lister[1]);
		iterator.all(basis1,new IteratorDirective(),lister[2]);
		assertEquals(list[1],list[2]);
		
		clearLists();
	}

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
