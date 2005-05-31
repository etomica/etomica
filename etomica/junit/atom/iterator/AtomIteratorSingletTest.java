package etomica.junit.atom.iterator;

import java.util.LinkedList;

import etomica.Atom;
import etomica.Space;
import etomica.atom.iterator.AtomIteratorSinglet;
import etomica.space2d.Space2D;


/**
 * Unit test for AtomIteratorSinglet class.
 *
 * @author David Kofke
 *
 */

/*
 * History
 * Created on May 30, 2005 by kofke
 */
public class AtomIteratorSingletTest extends IteratorTest {

    /**
     * 
     */
    public AtomIteratorSingletTest() {
        super();
    }
    
    public void setUp() {
        Space space = new Space2D();
        iterator = new AtomIteratorSinglet();
        testAtom1 = new Atom(space);
        testAtom2 = new Atom(space);
        list1 = makeTestList(new Atom[] {testAtom1});
        list2 = makeTestList(new Atom[] {testAtom2});
    }
    
    public void testIterator() {
        print("starting");
        LinkedList list = generalIteratorMethodTests(iterator);
        iterator.setAtom(testAtom1);
        list = generalIteratorMethodTests(iterator);
        assertEquals(list,list1);
        iterator.setAtom(null);
        assertNull(iterator.getAtom());
        list = generalIteratorMethodTests(iterator);
        assertNull(iterator.getAtom());
        assertTrue(list.size() == 0);
        iterator.setAtom(testAtom2);
        list = generalIteratorMethodTests(iterator);
        assertEquals(list, list2);
        iterator.setAtom(testAtom1);
        assertEquals(testAtom1, iterator.getAtom());
        list = generalIteratorMethodTests(iterator);
        assertEquals(list, list1);
        assertEquals(testAtom1, iterator.getAtom());
    }
    
    private AtomIteratorSinglet iterator;
    private Atom testAtom1, testAtom2;
    private LinkedList list1, list2;

}
