package etomica.junit.atom.iterator;

import java.util.LinkedList;

import etomica.atom.Atom;
import etomica.atom.AtomList;
import etomica.atom.iterator.ApiInterList;
import etomica.space.Space;
import etomica.space3d.Space3D;

/**
 * Unit test for ApiInterList
 * 
 * @author David Kofke
 *
 */

/*
 * History
 * Created on Jun 28, 2005 by kofke
 */
public class ApiInterListTest extends IteratorTestAbstract {

    public void testIterator() {
        
        ApiInterList api = new ApiInterList();
        Space space = Space3D.getInstance();
        
        //test that new iterator gives no iterates
        countTest(api, 0);
        
        AtomList outerList = new AtomList();
        AtomList innerList = new AtomList();
        api.setOuterList(outerList);
        api.setInnerList(innerList);
        
        //test both lists empty
        countTest(api, 0);
        
        //test inner list empty
        outerList.add(new Atom(space));
        countTest(api, 0);
        
        //test one atom in each list
        innerList.add(new Atom(space));
        countTest(api, 1);
        
        //populate outer list
        for(int i=2; i<5; i++) {
            outerList.add(new Atom(space));
            countTest(api, i);
        }
        
        //populate inner list
        for(int i=2; i<5; i++) {
            innerList.add(new Atom(space));
            countTest(api, 4*i);
        }

        //empty outer list
        outerList.clear();
        countTest(api, 0);
        //repopulate outer list
        for(int i=1; i<5; i++) {
            outerList.add(new Atom(space));
            countTest(api, 4*i);
        }
        
        //test constructor taking list arguments
        LinkedList list0 = countTest(api, 16);
        LinkedList list1 = countTest(new ApiInterList(outerList, innerList), 16);
        assertEquals(list0, list1);
        
        //test null lists
        api.setOuterList(null);
        countTest(api, 0);
        
        api.setOuterList(outerList);
        api.setInnerList(null);
        countTest(api, 0);
        
        api.setOuterList(null);
        countTest(api, 0);
        
        //test documented exceptions
        boolean exceptionThrown = false;
        try {
            new ApiInterList(outerList, outerList);
        } catch(IllegalArgumentException e) {exceptionThrown = true;}
        assertTrue(exceptionThrown);
        exceptionThrown = false;
        api.setOuterList(outerList);
        api.setInnerList(outerList);
        try {
            api.reset();
        } catch(IllegalStateException e) {exceptionThrown = true;}
        assertTrue(exceptionThrown);
    }
}
