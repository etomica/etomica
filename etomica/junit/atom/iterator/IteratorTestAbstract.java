package etomica.junit.atom.iterator;

import java.util.LinkedList;

import junit.framework.TestCase;
import etomica.atom.AtomPair;
import etomica.atom.AtomSet;
import etomica.atom.AtomsetArray;
import etomica.atom.IAtom;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomPairIterator;
import etomica.atom.iterator.AtomsetIterator;
import etomica.atom.iterator.IteratorDirective;
import etomica.junit.UnitTestUtil;

/**
 * Provides library methods to test the basic functioning of an iterator.
 * Most JUnit iterator tests extend this class.
 * 
 * @author Ken Benjamin and David Kofke
 *  
 */
public abstract class IteratorTestAbstract extends TestCase {


    /**
     * Clears the list of each lister in the array.
     */
    public void clearLists(Lister[] lister) {
        for (int i = 0; i < lister.length; i++) {
            lister[i].list.clear();
        }
    }

    /**
     * Prints all of the lists.
     */
    public void printLists(Lister[] lister) {
        if (UnitTestUtil.VERBOSE) {
            for (int i = 0; i < lister.length; i++) {
                System.out.println(lister[i]);
            }
            System.out.println();
        }
    }
    
    /**
     * Returns a list from the given atoms, which can be checked against
     * the list returned by generalIteratorMethodTests.
     */
    public LinkedList makeTestList(AtomSet[] atoms) {
        Lister lister = new Lister();
        for(int i=0; i<atoms.length; i++) {
            lister.actionPerformed(atoms[i]);
        }
        return lister.list;
    }
    
    protected void print(String string) {
        if(UnitTestUtil.VERBOSE) System.out.println(string);
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
    protected java.util.LinkedList generalIteratorMethodTests(
            AtomsetIterator iterator) {
        Lister[] lister = Lister.listerArray(5);

        //******* test of allAtoms
        iterator.allAtoms(lister[0]);
        iterator.allAtoms(lister[2]);

        printLists(lister);
        print("Just printed the lists at the beginning of generalIteratorMethodTests");
        assertEquals(lister[0].list, lister[2].list);

        //******* test of size
        print("Testing size -- iterator.size, lister.list.size: "+iterator.size()+" "+lister[0].list.size());
        assertEquals(iterator.size(), lister[0].list.size());

        AtomSet[] atoms = new AtomSet[lister[0].list.size()];
        
        //******* test of hasNext/next
        iterator.reset();
        while (iterator.hasNext()) {
            lister[3].actionPerformed(iterator.next());
        }
        assertEquals(lister[0].list, lister[3].list);

        //******* test of next
        iterator.reset();
        int j=0;
        while (iterator.hasNext()) {
            AtomSet nextAtom = iterator.next();
            if (nextAtom instanceof IAtom) {
                atoms[j] = nextAtom;
                if (!nextAtom.toString().equals(lister[0].list.get(j))) {
                    System.out.println(j+" "+nextAtom+" "+lister[0].list.get(j));
                }
                assertTrue(nextAtom.toString().equals(lister[0].list.get(j)));
            }
            else {
                atoms[j] = new AtomsetArray(nextAtom);
            }
            lister[1].actionPerformed(nextAtom);
            j++;
        }

        //******* test of nextAtom
        if(iterator instanceof AtomIterator) {
            iterator.reset();
            int i = 0;
            while(iterator.hasNext()) {
                IAtom next = ((AtomIterator)iterator).nextAtom();
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
            if (iterator != null) {
                assertNull(iterator.next());
            }
            assertFalse(iterator.hasNext());
        }
        print("Just tested unset method");

        //******* test of nBody
        iterator.reset();
        while (iterator.hasNext()) {
            assertEquals(iterator.next().count(), iterator.nBody());
        }

        print("Just tested nBody method");
        
//        if(iterator.nBody() == 2) {
//            iterator.reset();
//            while(iterator.hasNext()) {
//                AtomSet pair = iterator.next();
//                System.out.println(pair.getAtom(0).node.index() + " " + pair.getAtom(1).node.index());
//                assertTrue(pair.getAtom(0).node.index() < pair.getAtom(1).node.index());
//            }
//        }

        return lister[0].list;
    }

    /**
     * Tests that iterator gives no iterates
     */
    protected void testNoIterates(AtomsetIterator iterator) {
        LinkedList list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(), 0);
    }

    /**
     * Tests that iterates given by iterator match a specified list.
     * @param iterator the conditioned iterator
     * @param partners array of expected iterates
     * @return the Lister list of iterates
     */
    protected LinkedList testIterates(AtomIterator iterator, IAtom[] iterates) {
        LinkedList list = generalIteratorMethodTests(iterator);
        Lister test = new Lister();
        for(int i=0; i<iterates.length; i++) {
            test.actionPerformed(iterates[i]);
        }
        assertEquals(list, test.list);
        return list;
    }

    /**
     * Tests that iterates given by iterate match a specified list, for
     * iteration uplist or dnlist (not both)
     * @param iterator the conditioned iterate
     * @param iterate the atom0 expected in all pair iterates
     * @param partners array of atom1 expected in the pair iterates
     * @return the Lister list of iterates
     */
    protected LinkedList testApiIterates(AtomPairIterator iterator, IteratorDirective.Direction direction,
            IAtom iterate, IAtom[] partners) {
        LinkedList list = generalIteratorMethodTests(iterator);
        Lister test = new Lister();
        for(int i=0; i<partners.length; i++) {
            if(direction == IteratorDirective.Direction.UP) {
                test.actionPerformed(new AtomPair(iterate, partners[i]));
            } else {
                test.actionPerformed(new AtomPair(partners[i], iterate));
            }
        }
        assertEquals(list, test.list);
        return list;
    }
    
    /**
     * Same as testApiIterates, but with atom1 the same in all pair iterates, while
     * atom0 varies.
     */
    protected LinkedList testApiIteratesSwap(AtomPairIterator iterator, IAtom iterate, IAtom[] partners) {
        LinkedList list = generalIteratorMethodTests(iterator);
        Lister test = new Lister();
        for(int i=0; i<partners.length; i++) {
            test.actionPerformed(new AtomPair(partners[i],iterate));
        }
        assertEquals(list, test.list);
        return list;
    }


    /**
     * Tests that iterator gives two particular iterates
     */
    protected LinkedList testApiTwoIterates(AtomPairIterator iterator, AtomPair pair0, AtomPair pair1) {
        if((pair0.atom1 == null || pair0.atom0 == null) && (pair1.atom0 == null || pair1.atom1 == null)) {
            testNoIterates(iterator);
            return new LinkedList();
        }
        LinkedList list = generalIteratorMethodTests(iterator);
        Lister test = new Lister();
        test.actionPerformed(pair0);
        test.actionPerformed(pair1);
        assertEquals(list, test.list);
        return list;
    }

    /**
     * Tests that iterator gives a single particular iterate.
     */
    protected LinkedList testApiOneIterate(AtomPairIterator iterator, AtomPair pair) {
        if(pair.atom0 == null || pair.atom1 == null) {
            testNoIterates(iterator);
            return new LinkedList();
        }
        LinkedList list = generalIteratorMethodTests(iterator);
        Lister test = new Lister();
        test.actionPerformed(pair);
        assertEquals(list, test.list);
        return list;
    }

    /**
     * Tests that iterator gives a specific number of iterates.
     * @param api the iterator
     * @param n the number of iterates it should give
     * @return the Lister list of iterates
     */
    protected LinkedList countTest(AtomsetIterator api, int n) {
        LinkedList list = generalIteratorMethodTests(api);
        assertEquals(list.size(), n);
        return list;
    }

    /**
     * Tests that iterates given by iterate match a specified list, for
     * iteration uplist and then dnlist
     * @param iterator the conditioned iterate
     * @param iterate the atom0 expected in all pair iterates
     * @param up array of atom1 expected in the uplist pair iterates
     * @param dn array of atom1 expected in the dnlist pair iterates
     * @return the Lister list of iterates
     */
    protected LinkedList testApiIterates(AtomPairIterator iterator, IAtom iterate, IAtom[] up, IAtom[] dn) {
        LinkedList list = generalIteratorMethodTests(iterator);
        Lister test = new Lister();
        for(int i=0; i<up.length; i++) {
            test.actionPerformed(new AtomPair(iterate, up[i]));
        }
        for(int i=0; i<dn.length; i++) {
            test.actionPerformed(new AtomPair(dn[i], iterate));
        }
        assertEquals(list, test.list);
        return list;
    }

}