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
        Space space = Space2D.getInstance();
        singletIterator = new AtomIteratorSinglet();
        testAtom1 = new Atom(space);
        testAtom2 = new Atom(space);
        list1 = makeTestList(new Atom[] {testAtom1});
        list2 = makeTestList(new Atom[] {testAtom2});
    }
    
    public void testIterator() {
        print("starting");
        LinkedList list = generalIteratorMethodTests(singletIterator);
        singletIterator.setAtom(testAtom1);
        list = generalIteratorMethodTests(singletIterator);
        assertEquals(list,list1);
        singletIterator.setAtom(null);
        assertNull(singletIterator.getAtom());
        list = generalIteratorMethodTests(singletIterator);
        assertNull(singletIterator.getAtom());
        assertTrue(list.size() == 0);
        singletIterator.setAtom(testAtom2);
        list = generalIteratorMethodTests(singletIterator);
        assertEquals(list, list2);
        singletIterator.setAtom(testAtom1);
        assertEquals(testAtom1, singletIterator.getAtom());
        list = generalIteratorMethodTests(singletIterator);
        assertEquals(list, list1);
        assertEquals(testAtom1, singletIterator.getAtom());
    }
    
    private AtomIteratorSinglet singletIterator;
    private Atom testAtom1, testAtom2;
    private LinkedList list1, list2;

}
