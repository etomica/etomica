/*
 * Created on Oct 7, 2004
 */
package etomica.junit.atom.iterator;

import java.util.LinkedList;

import etomica.atom.AtomArrayList;
import etomica.atom.AtomLeaf;
import etomica.atom.iterator.AtomIteratorArrayListSimple;
import etomica.junit.UnitTestUtil;
import etomica.space3d.Space3D;

/**
 * Unit test for AtomIteratorArrayList.
 */
public class AtomIteratorArrayListTest extends IteratorTestAbstract {

    public AtomIteratorArrayListTest() {
        super();
        UnitTestUtil.VERBOSE = false;
    }
    
    public void testListVariations() {
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
            atomList.add(new AtomLeaf(Space3D.getInstance()));
        }
        list = generalIteratorMethodTests(iterator);
        
        //check that setList changes list
        AtomArrayList arrayList = new AtomArrayList();
        iterator.setList(arrayList);
        list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(), 0);
        arrayList.add(new AtomLeaf(Space3D.getInstance()));
        list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(), 1);
        
        //check handling of null list
        iterator.setList(null);
        list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(), 0);
        
    }

}

