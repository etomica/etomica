/*
 * Created on Oct 7, 2004
 */
package etomica.junit.atom.iterator;

import java.util.LinkedList;

import etomica.atom.Atom;
import etomica.atom.AtomArrayList;
import etomica.atom.iterator.AtomIteratorArrayListCompound;
import etomica.atom.iterator.AtomIteratorArrayListSimple;
import etomica.atom.iterator.AtomIteratorListSimple;
import etomica.junit.UnitTestUtil;
import etomica.space.Space;
import etomica.space1d.Space1D;

/**
 * Unit test for AtomIteratorListCompound.
 */
public class AtomIteratorListCompoundTest extends IteratorTestAbstract {

    public AtomIteratorListCompoundTest() {
        super();
        UnitTestUtil.VERBOSE = false;
    }
    
    public void testListVariations() {

        AtomIteratorArrayListCompound iterator = new AtomIteratorArrayListCompound();
        
        //make sure new iterator gives no iterates
        LinkedList list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(), 0);
        // make empty lists to start
        AtomArrayList atomList1 = new AtomArrayList();
        AtomArrayList atomList2 = new AtomArrayList();
        AtomArrayList atomList3 = new AtomArrayList();
        AtomArrayList[] array = new AtomArrayList[] {atomList1, atomList2, atomList3};
        iterator.setLists(array);
        
        //test various combinations of lists of different sizes
        doVariations(iterator, array, 2);
        
        atomList1.remove(0);
        atomList1.remove(0);
        atomList2.remove(0);
        
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
        iterator.setLists(new AtomArrayList[] {atomList1});
        list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(),i1);
        
        //add same list twice
        iterator.setLists(new AtomArrayList[] {atomList1, atomList1});
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
    
    public void doVariations(AtomIteratorArrayListCompound iterator, AtomArrayList[] lists, int varIndex) {
        
        Space space = Space1D.getInstance();
        
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
    
    private Atom[] makeArray(AtomArrayList[] lists) {
        AtomArrayList aList = new AtomArrayList();
        AtomIteratorArrayListSimple iter = new AtomIteratorArrayListSimple();
        for(int i=0; i<lists.length; i++) {
            iter.setList(lists[i]);
            iter.reset();
            while (iter.hasNext()) {
                aList.add(iter.nextAtom());
            }
        }
        return aList.toArray();
    }

}

