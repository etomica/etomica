/*
 * Created on Oct 7, 2004
 */
package etomica.junit.atom.iterator;

import java.util.LinkedList;

import etomica.atom.Atom;
import etomica.atom.AtomArrayList;
import etomica.atom.iterator.AtomIteratorArrayListSimple;
import etomica.junit.UnitTest;
import etomica.space.Space;
import etomica.space1d.Space1D;

/**
 * Unit test for AtomIteratorArrayList.
 */
public class AtomIteratorArrayListTest extends IteratorTest {

    public AtomIteratorArrayListTest() {
        super();
        UnitTest.VERBOSE = false;
    }
    
    public void testListVariations() {
        Space space = Space1D.getInstance();
        AtomIteratorArrayListSimple iterator = new AtomIteratorArrayListSimple();
        
        //make sure new iterator gives no iterates
        LinkedList list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(), 0);
        
        // make empty list to start
        AtomArrayList atomList = new AtomArrayList();
        iterator.setList(atomList);
        
        //add some atoms and check each time
        for(int i=0; i<=10; i++) {
            list = generalIteratorMethodTests(iterator);
            assertEquals(list.size(), i);
            atomList.add(new Atom(space));
        }
        list = generalIteratorMethodTests(iterator);
        
        //check that setList changes list
        AtomArrayList arrayList = new AtomArrayList();
        iterator.setList(arrayList);
        list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(), 0);
        arrayList.add(new Atom(space));
        list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(), 1);
        
        //check handling of null list
        iterator.setList(null);
        list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(), 0);
        
    }

}

