package etomica.junit;

import java.util.Collections;
import java.util.LinkedList;

import etomica.*;
import etomica.atom.AtomList;
import etomica.atom.iterator.AtomIteratorList;
import junit.framework.*;

/**
 * @author KMB
 *
 */
public abstract class IteratorTest extends TestCase {

	/**
	 * Clears all of the lists.
	*/ 
		public void clearLists(Lister[] lister) {
			for (int i=0;i<lister.length;i++) {
				lister[i].list.clear();
			}
		}
		
	/**
	 * Prints all of the lists.
	 */
		public void printLists(Lister[] lister) {
			for (int i=0;i<lister.length;i++) {
				System.out.println(lister[i].list);
			}
			System.out.println();
		}

/**
 * this method tests the different methods associated with a list iterator,
 * such as hasNext, next, contains, size, etc.
 * @param iterator
 * @return
 */	
	public java.util.LinkedList generalIteratorMethodTests(AtomIterator iterator) {
		// initialize lists here with code or separate setUp method
		Lister[] lister = Lister.listerArray(4);
		iterator.allAtoms(lister[0]);
		iterator.allAtoms(lister[2]);
		
		printLists(lister);
		System.out.println("Just printed the lists at the beginning of generalIteratorMethodTests");
// 		test whether iterator does same thing twice		
		assertEquals(lister[0].list, lister[2].list);
//		clearLists(lister);

		iterator.reset();
		while(iterator.hasNext()) {
			AtomSet peekAtom = iterator.peek();
			assertEquals(peekAtom, iterator.next());
			lister[1].actionPerformed(peekAtom);
		}
		
		iterator.reset();
		while(iterator.hasNext()) {
			lister[3].actionPerformed(iterator.nextAtom());
		}
		assertEquals(lister[0].list, lister[3].list);
		
		//test that allAtoms and hasNext/next give same set of iterates
		assertEquals(lister[0].list, lister[1].list);
		System.out.println("Just tested for allAtoms and hasNext/next");
		printLists(lister);
//		clearLists(lister);
		
		//test operation of unset method
		iterator.reset();
		iterator.unset();
		assertFalse(iterator.hasNext());
//		assertNull(iterator.next()[0]);
		assertNull(iterator.next());
		assertFalse(iterator.hasNext());
		for(int i=0; i<5; i++) {
			if(iterator != null) assertNull(iterator.nextAtom());
			assertFalse(iterator.hasNext());
		}
		System.out.println("Just tested unset method");
		
		//test size method
		assertEquals(iterator.size(), lister[0].list.size());
		System.out.println("Just tested size method");
		
		//test contains method
		if(iterator != null) {
			AtomList atomList = new AtomList(iterator);
			AtomIteratorList listIterator = new AtomIteratorList(atomList);
			listIterator.reset();
			while(listIterator.hasNext()) {
				AtomSet nextAtom=listIterator.next();
				System.out.println("nextAtom = "+ nextAtom);
				assertTrue(iterator.contains(nextAtom));
				System.out.println("Contains works with hasNext method, contains equals "+ iterator.contains(nextAtom));
			}
			assertFalse(iterator.contains(null));
//			System.out.println("Contains equals "+ iterator.contains(nextAtom));

		}
		System.out.println("Just tested contains method");
		
		//test nBody
		iterator.reset();

//		assertEquals(iterator.next().length, iterator.nBody());

// 		Error occurred here, kmb 4/19/05		
// 		Added if iterator.hasNext conditional to handle null pointer problems when
// 		the iterator was empty.  kmb 4/27/05
		if (iterator.hasNext()) {
			assertEquals(iterator.next().count(), iterator.nBody());
		}
		
		System.out.println("Just tested nBody method");
		
		return lister[0].list;
	}

	
}
