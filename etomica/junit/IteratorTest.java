package etomica.junit;

import junit.framework.TestCase;
import etomica.Atom;
import etomica.AtomIterator;
import etomica.AtomSet;
import etomica.AtomsetIterator;

/**
 * Provides library methods to test the basic functioning of an iterator.
 * 
 * @author KMB
 *  
 */
public class IteratorTest extends TestCase {

    /**
     * Declare constructor private to prevent instantiation.
     */
    private IteratorTest() {
    }

    /**
     * Clears the list of each lister in the array.
     */
    public static void clearLists(Lister[] lister) {
        for (int i = 0; i < lister.length; i++) {
            lister[i].list.clear();
        }
    }

    /**
     * Prints all of the lists.
     */
    public static void printLists(Lister[] lister) {
        if (UnitTest.VERBOSE) {
            for (int i = 0; i < lister.length; i++) {
                System.out.println(lister[i].list);
            }
            System.out.println();
        }
    }
    
    private static void print(String string) {
        if(UnitTest.VERBOSE) System.out.println(string);
    }

    /**
     * Tests all the elementary methods of an atomset iterator. Methods tested are
     * allAtoms, hasNext/next, peek, size, contains, nBody, unset, reset.  Suitable
     * for use by any AtomsetIterator.  Tests nextAtom if iterator is instance
     * of AtomIterator
     * 
     * @param iterator
     *            an iterator in a condition to be tested; reset and iteration
     *            will be performed repeatedly on the iterator
     * @return a list of Lister instances for the atoms given by the iterator
     */
    public static java.util.LinkedList generalIteratorMethodTests(
            AtomsetIterator iterator) {
        Lister[] lister = Lister.listerArray(5);

        //******* test of allAtoms
        iterator.allAtoms(lister[0]);
        iterator.allAtoms(lister[2]);

        printLists(lister);
        print("Just printed the lists at the beginning of generalIteratorMethodTests");
        assertEquals(lister[0].list, lister[2].list);

        //******* test of size
        assertEquals(iterator.size(), lister[0].list.size());
        print("Just tested size method");

        AtomSet[] atoms = new AtomSet[lister[0].list.size()];
        
        //******* test of peek
        iterator.reset();
        int j=0;
        while (iterator.hasNext()) {
            AtomSet peekAtom = iterator.peek();
            atoms[j++] = peekAtom;
            assertEquals(peekAtom, iterator.next());
            lister[1].actionPerformed(peekAtom);
        }

        //******* test of hasNext/next
        iterator.reset();
        while (iterator.hasNext()) {
            lister[3].actionPerformed(iterator.next());
        }
        assertEquals(lister[0].list, lister[3].list);

        //******* test of nextAtom
        if(iterator instanceof AtomIterator) {
            iterator.reset();
            int i = 0;
            while(iterator.hasNext()) {
                Atom next = ((AtomIterator)iterator).nextAtom();
                assertEquals(next, atoms[i++]);
            }
        }
        
        //******* test that allAtoms and hasNext/next give same set of iterates
        assertEquals(lister[0].list, lister[1].list);
        print("Just tested for allAtoms and hasNext/next");
        printLists(lister);

        //******* test of unset
        iterator.reset();
        iterator.unset();
        assertFalse(iterator.hasNext());
        assertNull(iterator.next());
        assertFalse(iterator.hasNext());
        //******* test that next calls cannot cause hasNext to become true
        //******* test that iterate is null for hasNext false
        for (int i = 0; i < 5; i++) {
            if (iterator != null)
                assertNull(iterator.next());
            assertFalse(iterator.hasNext());
        }
        print("Just tested unset method");


        //******* test of contains
        if (iterator != null) {
            for(int i=0; i<atoms.length; i++) {
                AtomSet nextAtom = atoms[i];
                print("nextAtom = " + nextAtom);
                assertTrue(iterator.contains(nextAtom));
                print("Contains works with hasNext method, contains equals "
                                    + iterator.contains(nextAtom));
            }
            assertFalse(iterator.contains(null));
        }
        print("Just tested contains method");

        //******* test of nBody
        iterator.reset();
        if (iterator.hasNext()) {
            assertEquals(iterator.next().count(), iterator.nBody());
        }

        print("Just tested nBody method");

        return lister[0].list;
    }

}