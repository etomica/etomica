package etomica.junit;

import etomica.*;
import etomica.action.AtomActionAdapter;

import java.util.LinkedList;

/**
 * @author Andrew Walker
 * 
 * This class test the methods of the AtomIteratorListSimple class using JUnit.
 */

public class AtomIteratorListSimpleTest extends IteratorTest {

	private AtomIteratorListSimple iterator;
	private LinkedList list3;

	public static void main(String[] args) {
		junit.textui.TestRunner.run(AtomIteratorListSimpleTest.class);
	}

/**
 * This method is called before any test is run and sets up the basics of the
 * etomica simulation and testing elements.
 */
	protected void setUp() {
		super.setUp();
		iterator = new AtomIteratorListSimple(atomList1);
	}

/**
 * Tests the size() method of the iterator.
 */
	public void testSize() {
		assertEquals(iterator.size(),nMolecules);
	}
	
	public void testUnSet() {
		
		iterator.reset();
		iterator.unset();
		assertEquals(false,iterator.hasNext());
	}
	
	public void testReset() {
		
		iterator.reset();
		Atom a = iterator.next();
		iterator.next();
		iterator.reset();
		assertEquals(a,iterator.next());
	}
	
	public void testNextLoop() {
		
		iterator.reset();
		
		for (int i=0;i<atomList1.size();i++) {
			lister[0].actionPerformed(atomList1.get(i));
		}
		while (iterator.hasNext()) {
			lister[1].actionPerformed(iterator.next());
		}
		assertEquals(lister[0].list,lister[1].list);

		clearLists();
	}
	
	public void testAllMethods() {
		
		iterator.reset();
		IteratorDirective dummy = new IteratorDirective();
		list3 = new LinkedList();
		AtomActionAdapter lister3 = new AtomActionAdapter() {
			public void actionPerformed(Atom a) {
				list3.add(a.toString());
			}
		};

		for (int i=0;i<atomList1.size();i++) {
			lister[0].actionPerformed(atomList1.get(i));
		}
		
		iterator.all(basis1,dummy,lister[1]);
		iterator.all(atomList1,dummy,lister[2]);
		iterator.allAtoms(lister3);
		
		assertEquals(lister[0].list,lister[1].list);
		assertEquals(lister[1].list,lister[2].list);
		assertEquals(lister[2].list,list3);
		
		clearLists();
	}
	
	public void testPeek() {
		
		iterator.reset();
		
		assertEquals(first1,iterator.peek());
		assertEquals(first1,iterator.peek());
	}
	
	public void testContains() {
		
		iterator.reset();
		
		assertEquals(true,iterator.contains(first1));
		assertEquals(false,iterator.contains(outsideAtom));
	}
	
}