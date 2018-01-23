/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/*
 * Created on Oct 7, 2004
 */
package etomica.atom.iterator;

import java.util.LinkedList;

import etomica.atom.Atom;
import etomica.atom.AtomArrayList;
import etomica.UnitTestUtil;
import etomica.space3d.Space3D;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import static etomica.atom.iterator.IteratorTestAbstract.generalIteratorMethodTests;

/**
 * Unit test for AtomIteratorArrayList.
 */
class AtomIteratorArrayListTest {

    public AtomIteratorArrayListTest() {
        super();
        UnitTestUtil.VERBOSE = false;
    }
    
    @Test
    public void testListVariations() {
        AtomIteratorArrayListSimple iterator = new AtomIteratorArrayListSimple();
        
        //make sure new iterator gives no iterates
        LinkedList list = generalIteratorMethodTests(iterator);
        Assertions.assertEquals(list.size(), 0);
        
        // make empty list to start
        AtomArrayList atomList = new AtomArrayList();
        iterator.setList(atomList);
        
        //add some atoms and check each time
        for(int i=0; i<=10; i++) {
            list = generalIteratorMethodTests(iterator);
            Assertions.assertEquals(list.size(), i);
            atomList.add(new Atom(Space3D.getInstance()));
        }
        list = generalIteratorMethodTests(iterator);
        
        //check that setList changes list
        AtomArrayList arrayList = new AtomArrayList();
        iterator.setList(arrayList);
        list = generalIteratorMethodTests(iterator);
        Assertions.assertEquals(list.size(), 0);
        arrayList.add(new Atom(Space3D.getInstance()));
        list = generalIteratorMethodTests(iterator);
        Assertions.assertEquals(list.size(), 1);
        
        //check handling of null list
        iterator.setList(null);
        boolean exceptionThrown = false;
        try {
            list = generalIteratorMethodTests(iterator);
        }
        catch (RuntimeException e) {
            exceptionThrown = true;
        }
        Assertions.assertTrue(exceptionThrown);
        
    }

}

