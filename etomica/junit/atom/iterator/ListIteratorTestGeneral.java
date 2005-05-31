/*
 * Created on Oct 1, 2004
 */
package etomica.junit.atom.iterator;
import java.util.Collections;
import java.util.LinkedList;

import etomica.Atom;
import etomica.IteratorDirective;
import etomica.atom.AtomLinker;
import etomica.atom.AtomList;
import etomica.atom.iterator.AtomIteratorList;
import etomica.atom.iterator.AtomsetIteratorListDependent;
import etomica.junit.UnitTest;

/**
 * Tests a general list iterator, ensuring that effect of direction, first,
 * and terminator are correct, and checking for appropriate behavior with 
 * or without tabs in list.
 * 
 *  @author Ken Benjamin
 *
 */
public class ListIteratorTestGeneral extends ListIteratorTest {
		
	protected final int nLists=10;
	protected LinkedList[] lists = new LinkedList[nLists];
	protected AtomList atomListTemp;
	protected LinkedList atomListUpFirst, atomListDownLast;
	protected AtomLinker firstLinker,lastLinker,midLinker;
	protected Atom firstAtom,middleAtom,lastAtom;


	public ListIteratorTestGeneral() {
		super(new AtomIteratorList());
	}
	
	public void setUp() {
		//does nothing
	}
	
		
	public void iteratorStateTests(AtomsetIteratorListDependent listIterator) {
		AtomIteratorList iterator = (AtomIteratorList)listIterator;
		iterator.setTerminatorType(AtomLinker.Tab.HEADER_TAB);
		listElements(listIterator.getList());
	 	
		/**  Sets the iteration direction to UP and test the methods in turn with
		 * 	 the first element equal to first, middle, and last.
		 */
		iterator.setDirection(IteratorDirective.UP);
		if(UnitTest.VERBOSE) System.out.println("the current iterator direction is: "+ iterator.getDirection());
		assertEquals(iterator.getDirection(), IteratorDirective.UP);
		// set the first element in the iterator to be the first atom in the list
		iterator.setFirst(firstAtom);
		if(UnitTest.VERBOSE) System.out.println("iterator.setFirst equals first");
		lists[0] = generalIteratorMethodTests(iterator);
		atomListUpFirst = lists[0];
		checkFirstLast(lists[0], firstAtom, lastAtom);
		
		AtomIteratorList copyIterator = new AtomIteratorList(iterator);
		LinkedList copyList = generalIteratorMethodTests(copyIterator);
		assertEquals(lists[0], copyList);
//		iteratorUpFirst=iterator;
//		iterator.allAtoms(lister[5]);
//		printLists(lister);
		// set the first element in the iterator to be the middle atom in the list
		iterator.setFirst(middleAtom);
		if(UnitTest.VERBOSE) System.out.println("iterator.setFirst equals middle");
		lists[1] = generalIteratorMethodTests(iterator);
		checkFirstLast(lists[1], middleAtom, lastAtom);
//		printLists(lister);
		// set the first element in the iterator to be the last atom in the list
		iterator.setFirst(lastAtom);
		if(UnitTest.VERBOSE) System.out.println("iterator.setFirst equals last");
		lists[2] = generalIteratorMethodTests(iterator);
		checkFirstLast(lists[2], lastAtom, lastAtom);
		///		printLists(lister);
		// set the first element in the iterator with no argument to verify it's
		// default behavior is the same as with setFirst equal to the first atom
		// in the list
/*		iterator.setFirst();  // won't accept setFirst with no argument
		System.out.println("iterator.setFirst equals last");
		generalIteratorMethodTests(iterator);
		iteratorUpFirstNoArg=iterator;
*///		printLists(lister);
		
		//Setting new list should keep direction same, but make header the first and last
		iterator.setList(iterator.getList());
		lists[3] = generalIteratorMethodTests(iterator);
		assertEquals(lists[0], lists[3]);
		
		iterator.setFirst(lastAtom);
		iterator.setFirst((AtomLinker)null);
		lists[3] = generalIteratorMethodTests(iterator);
		assertEquals(lists[0], lists[3]);
		
		/**  Sets the iteration direction to DOWN and test the methods in turn with
		 * 	 the first element equal to first, middle, and last.
		 */
		iterator.setDirection(IteratorDirective.DOWN);
		assertEquals(IteratorDirective.DOWN, iterator.getDirection());
		if(UnitTest.VERBOSE) System.out.println("the current iterator direction is: "+ iterator.getDirection());
		// set the first element in the iterator to be the first atom in the list
		iterator.setFirst(firstAtom);
		if(UnitTest.VERBOSE) System.out.println("iterator.setFirst equals first");
		lists[3] = generalIteratorMethodTests(iterator);
		checkFirstLast(lists[3], firstAtom, firstAtom);
//		iteratorDownFirst=iterator;
//		iterator.allAtoms(lister[6]);
//		printLists(lister);
		// set the first element in the iterator to be the middle atom in the list
		iterator.setFirst(middleAtom);
		if(UnitTest.VERBOSE) System.out.println("iterator.setFirst equals middle");
		lists[4] = generalIteratorMethodTests(iterator);
		checkFirstLast(lists[4], middleAtom, firstAtom);

//		printLists(lister);
		// set the first element in the iterator to be the last atom in the list
		iterator.setFirst(lastAtom);
		if(UnitTest.VERBOSE) System.out.println("iterator.setFirst equals last");
		lists[5] = generalIteratorMethodTests(iterator);
		checkFirstLast(lists[5],lastAtom, firstAtom);

//		iteratorDownLast=iterator;
//		printLists(lister);
//		iterator.allAtoms(lister[3]);
		atomListDownLast=lists[5];
//		atomListDownLast=lister[0].list;
		Collections.reverse(atomListDownLast);
		assertEquals(atomListUpFirst, atomListDownLast);
		if(UnitTest.VERBOSE) System.out.println("Just tested with reverse collections and compared listUpFirst and listDownLast");
//		printLists(lister);
//		clearLists(lister);
//		Commented out all calls to printLists and clearLists due to error with lister, 12/06/04	
		iterator.setDirection(IteratorDirective.UP);
		iterator.setTerminatorType(tabType1);
		lists[6] = generalIteratorMethodTests(iterator);
		iterator.setTerminatorType(tabType2);
		lists[7] = generalIteratorMethodTests(iterator);
		
	}
	
	private void checkFirstLast(LinkedList list, Atom first, Atom last) {
		if(list.size() > 0) {
			assertEquals(list.getFirst(),first.toString());
			assertEquals(list.getLast(), last.toString());
		}
	}
		

	
	public void listElements(AtomList atomList) {
		
		firstAtom = atomList.getFirst();
		lastAtom = atomList.getLast();
		int mid=atomList.size()/2;
		// check list size > 0, else get first one
		middleAtom = atomList.size() > 0 ? atomList.get(mid) : firstAtom;
		firstLinker = atomList.firstEntry();
		lastLinker = atomList.lastEntry();
		midLinker = atomList.size() > 0 ? atomList.entry(mid) : firstLinker;	
	}



}
		
