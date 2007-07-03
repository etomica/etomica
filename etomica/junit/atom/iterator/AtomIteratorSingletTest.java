package etomica.junit.atom.iterator;

import java.util.LinkedList;

import etomica.atom.AtomLeaf;
import etomica.atom.AtomSet;
import etomica.atom.AtomSetSinglet;
import etomica.atom.IAtom;
import etomica.atom.iterator.AtomIteratorSinglet;
import etomica.space3d.Space3D;


/**
 * Unit test for AtomIteratorSinglet class.
 *
 * @author David Kofke
 *
 */
public class AtomIteratorSingletTest extends IteratorTestAbstract {

    public AtomIteratorSingletTest() {
        super();
    }
    
    public void setUp() {
        singletIterator = new AtomIteratorSinglet();
        testAtom1 = new AtomLeaf(Space3D.getInstance());
        testAtom2 = new AtomLeaf(Space3D.getInstance());
        list1 = makeTestList(new AtomSet[] {new AtomSetSinglet(testAtom1)});
        list2 = makeTestList(new AtomSet[] {new AtomSetSinglet(testAtom2)});
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
    private IAtom testAtom1, testAtom2;
    private LinkedList list1, list2;

}
