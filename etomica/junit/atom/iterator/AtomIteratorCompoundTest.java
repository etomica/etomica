/*
 * Created on Oct 7, 2004
 */
package etomica.junit.atom.iterator;

import java.util.LinkedList;

import etomica.atom.Atom;
import etomica.atom.AtomList;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorCompound;
import etomica.atom.iterator.AtomIteratorListSimple;
import etomica.junit.UnitTestUtil;
import etomica.space.Space;
import etomica.space1d.Space1D;

/**
 * Unit test for AtomIteratorListCompound.
 */
public class AtomIteratorCompoundTest extends IteratorTestAbstract {

    public AtomIteratorCompoundTest() {
        super();
        UnitTestUtil.VERBOSE = false;
    }
    
    public void testListVariations() {

        AtomIteratorCompound iterator = new AtomIteratorCompound();
        
        //make sure new iterator gives no iterates
        LinkedList list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(), 0);
        // make empty lists to start
        AtomList atomList1 = new AtomList();
        AtomList atomList2 = new AtomList();
        AtomList atomList3 = new AtomList();
        AtomIteratorListSimple atomIterator1 = new AtomIteratorListSimple(atomList1);
        AtomIteratorListSimple atomIterator2 = new AtomIteratorListSimple(atomList2);
        AtomIteratorListSimple atomIterator3 = new AtomIteratorListSimple(atomList3);
        
        AtomIteratorListSimple[] array = new AtomIteratorListSimple[] {atomIterator1, atomIterator2, atomIterator3};
        iterator.setIterators(array);
        
        //test various combinations of lists of different sizes
        doVariations(iterator, array, 2);
        
        atomList1.removeFirst();
        atomList1.removeFirst();
        atomList2.removeFirst();
        
        int i1 = atomList1.size();
        int i2 = atomList2.size();
        int i3 = atomList3.size();

        //remove list
        iterator.removeIterator(atomIterator1);
        list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(),i2+i3);
        
        //replace list
        iterator.addIterator(atomIterator1);
        list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(),i1+i2+i3);

        //add duplicate of list
        iterator.addIterator(atomIterator2);
        list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(),i1+2*i2+i3);

        //setIterators
        iterator.setIterators(new AtomIterator[] {atomIterator1});
        list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(),i1);
        
        //add same list twice
        iterator.setIterators(new AtomIterator[] {atomIterator1, atomIterator1});
        list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(),2*i1);
        
        //remove duplicate reference to list
        iterator.removeIterator(atomIterator1);
        list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(),i1);
        
        //add null list
        iterator.addIterator(null);
        list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(),i1);
        
        //remove list that wasn't added
        iterator.removeIterator(atomIterator2);
        list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(),i1);
        
        iterator.setIterators(null);
        list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(), 0);
        
    }
    
    public void doVariations(AtomIteratorCompound iterator, AtomIteratorListSimple[] iterators, int varIndex) {
        
        Space space = Space1D.getInstance();
        
        AtomList[] lists = new AtomList[iterators.length];
        for(int i=0; i<lists.length; i++) {
            lists[i] = iterators[i].getList();
        }
        //add some atoms and check each time
        lists[varIndex].clear();
        for(int i=0; i<=4; i++) {
            if(varIndex > 0) {
                doVariations(iterator, iterators, varIndex-1);
            } else {
                testIterates(iterator, makeArray(lists));
            }
            lists[varIndex].add(new Atom(space));
        }
        generalIteratorMethodTests(iterator);
    }
    
    private Atom[] makeArray(AtomList[] lists) {
        AtomList aList = new AtomList();
        AtomIteratorListSimple iter = new AtomIteratorListSimple();
        for(int i=0; i<lists.length; i++) {
            iter.setList(lists[i]);
            aList.addAll(iter);
        }
        return aList.toArray();
    }

}

