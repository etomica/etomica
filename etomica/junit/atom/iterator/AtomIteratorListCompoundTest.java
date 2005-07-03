/*
 * Created on Oct 7, 2004
 */
package etomica.junit.atom.iterator;

import java.util.LinkedList;

import etomica.Atom;
import etomica.Space;
import etomica.atom.AtomList;
import etomica.atom.iterator.AtomIteratorListCompound;
import etomica.atom.iterator.AtomIteratorListSimple;
import etomica.junit.UnitTest;
import etomica.space1d.Space1D;

/**
 * Unit test for AtomIteratorListCompound.
 */
public class AtomIteratorListCompoundTest extends IteratorTest {

    public AtomIteratorListCompoundTest() {
        super();
        UnitTest.VERBOSE = false;
    }
    
    public void testListVariations() {

        AtomIteratorListCompound iterator = new AtomIteratorListCompound();
        
        //make sure new iterator gives no iterates
        LinkedList list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(), 0);
        // make empty lists to start
        AtomList atomList1 = new AtomList();
        AtomList atomList2 = new AtomList();
        AtomList atomList3 = new AtomList();
        AtomList[] array = new AtomList[] {atomList1, atomList2, atomList3};
        iterator.setLists(array);
        
        //test various combinations of lists of different sizes
        doVariations(iterator, array, 2);
        
        atomList1.removeFirst();
        atomList1.removeFirst();
        atomList2.removeFirst();
        
        int i1 = atomList1.size();
        int i2 = atomList2.size();
        int i3 = atomList3.size();

        //remove list
        iterator.removeList(atomList1);
        list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(),i2+i3);
        
        //replace list
        iterator.addList(atomList1);
        list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(),i1+i2+i3);

        //add duplicate of list
        iterator.addList(atomList2);
        list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(),i1+2*i2+i3);

        //setLists
        iterator.setLists(new AtomList[] {atomList1});
        list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(),i1);
        
        //add same list twice
        iterator.setLists(new AtomList[] {atomList1, atomList1});
        list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(),2*i1);
        
        //remove duplicate reference to list
        iterator.removeList(atomList1);
        list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(),i1);
        
        //add null list
        iterator.addList(null);
        list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(),i1);
        
        //remove list that wasn't added
        iterator.removeList(atomList2);
        list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(),i1);
        
        iterator.setLists(null);
        list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(), 0);
        
    }
    
    public void doVariations(AtomIteratorListCompound iterator, AtomList[] lists, int varIndex) {
        
        Space space = new Space1D();
        
        //add some atoms and check each time
        lists[varIndex].clear();
        for(int i=0; i<=4; i++) {
            if(varIndex > 0) {
                doVariations(iterator, lists, varIndex-1);
            } else {
//                for(int k=0; k<lists.length; k++) System.out.print(lists[k].size()+" ");
//                System.out.println();
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

