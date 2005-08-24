package etomica.junit.atom.iterator;

import java.util.LinkedList;

import etomica.Space;
import etomica.atom.Atom;
import etomica.atom.AtomPair;
import etomica.atom.AtomSet;
import etomica.atom.iterator.ApiSinglet;
import etomica.space2d.Space2D;

/**
 * Unit test for ApiSinglet class.
 * 
 * @author David Kofke
 *  
 */

/*
 * History Created on May 30, 2005 by kofke
 */
public class ApiSingletTest extends IteratorTest {

    /**
     *  
     */
    public ApiSingletTest() {
        super();
    }

    public void setUp() {
        Space space = Space2D.getInstance();
        iterator = new ApiSinglet();
        testAtom1 = new Atom(space);
        testAtom2 = new Atom(space);
        testAtom3 = new Atom(space);
        list1 = makeTestList(new AtomSet[] { new AtomPair(testAtom1, testAtom2) });
        list2 = makeTestList(new AtomSet[] { new AtomPair(testAtom2, testAtom1) });
        list3 = makeTestList(new AtomSet[] { new AtomPair(testAtom2, testAtom3) });
    }

    public void testIterator() {

        LinkedList list = generalIteratorMethodTests(iterator);
        assertTrue(list.size() == 0);

        iterator.setPair(testAtom1, testAtom2);
        list = generalIteratorMethodTests(iterator);
        assertEquals(list, list1);
        assertFalse(list.equals(list2));
        assertEquals(list.size(), 1);

        iterator.setPair(null, null);
        assertNotNull(iterator.getPair());
        list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(), 0);

        iterator.setPair(null, testAtom2);
        assertNotNull(iterator.getPair());
        list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(), 0);

        iterator.setPair(testAtom1, null);
        assertNotNull(iterator.getPair());
        list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(), 0);

        iterator.setPair(testAtom2, testAtom1);
        assertFalse(new AtomPair(testAtom1, testAtom2).equals(iterator.getPair()));
        assertEquals(new AtomPair(testAtom2, testAtom1), iterator.getPair());
        assertFalse(iterator.contains(new AtomPair(testAtom1, testAtom2)));
        assertFalse(iterator.contains(new AtomPair(testAtom2, testAtom3)));
        assertTrue(iterator.contains(new AtomPair(testAtom2, testAtom1)));

        iterator.getPair().atom1 = null;
        list = generalIteratorMethodTests(iterator);
        assertFalse(list.equals(list1));
        
        iterator.getPair().atom1 = testAtom3;
        list = generalIteratorMethodTests(iterator);
        assertEquals(list, list3);
    }

    private ApiSinglet iterator;
    private Atom testAtom1, testAtom2, testAtom3;
    private LinkedList list1, list2, list3;

}