/*
 * Created on Oct 7, 2004
 */
package etomica.junit.atom.iterator;

import java.util.LinkedList;

import etomica.atom.Atom;
import etomica.atom.AtomList;
import etomica.atom.iterator.AtomIteratorListSimple;
import etomica.junit.UnitTestUtil;
import etomica.space.Space;
import etomica.space1d.Space1D;

/**
 * Unit test for AtomIteratorListSimple.
 */
public class AtomIteratorListSimpleTest extends IteratorTestAbstract {

    public AtomIteratorListSimpleTest() {
        super();
        UnitTestUtil.VERBOSE = false;
    }
    
    public void testListVariations() {
        Space space = Space1D.getInstance();
        AtomIteratorListSimple iterator = new AtomIteratorListSimple();
        
        //make sure new iterator gives no iterates
        LinkedList list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(), 0);
        
        // make empty list to start
        AtomList atomList = new AtomList();
        iterator.setList(atomList);
        assertEquals(atomList, iterator.getList());
        
        //add some atoms and check each time
        for(int i=0; i<=10; i++) {
            list = generalIteratorMethodTests(iterator);
            assertEquals(list.size(), i);
            atomList.add(new Atom());
        }
        list = generalIteratorMethodTests(iterator);
        
        //check that iterator functions ok if atoms are deleted
        //as they are being iterated via hasNext/next
        int n = atomList.size();
        Lister lister = new Lister();
        iterator.reset();
        while(iterator.hasNext()) {
            Atom deletedAtom = iterator.nextAtom();
            atomList.remove(deletedAtom);
            n--;
            lister.actionPerformed(deletedAtom);
        }
        assertEquals(n, 0);
        assertEquals(list, lister.list);
        
        //check that setList changes list
        iterator.setList(new AtomList());
        list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(), 0);
        iterator.getList().add(new Atom());
        list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(), 1);
        
        iterator.setList(null);
        list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(), 0);
        
    }

}

