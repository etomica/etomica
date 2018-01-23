/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom.iterator;

import java.util.LinkedList;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.atom.Atom;
import etomica.atom.AtomSetSinglet;
import etomica.space3d.Space3D;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import static etomica.atom.iterator.IteratorTestAbstract.*;


/**
 * Unit test for AtomIteratorSinglet class.
 *
 * @author David Kofke
 *
 */
class AtomIteratorSingletTest {

    public AtomIteratorSingletTest() {
        super();
    }
    
    @BeforeEach
    public void setUp() {
        singletIterator = new AtomIteratorSinglet();
        testAtom1 = new Atom(Space3D.getInstance());
        testAtom2 = new Atom(Space3D.getInstance());
        list1 = makeTestList(new IAtomList[] {new AtomSetSinglet(testAtom1)});
        list2 = makeTestList(new IAtomList[] {new AtomSetSinglet(testAtom2)});
    }
    
    @Test
    public void testIterator() {
        print("starting");
        LinkedList list = generalIteratorMethodTests(singletIterator);
        singletIterator.setAtom(testAtom1);
        list = generalIteratorMethodTests(singletIterator);
        Assertions.assertEquals(list, list1);
        singletIterator.setAtom(null);
        Assertions.assertNull(singletIterator.getAtom());
        list = generalIteratorMethodTests(singletIterator);
        Assertions.assertNull(singletIterator.getAtom());
        Assertions.assertTrue(list.size() == 0);
        singletIterator.setAtom(testAtom2);
        list = generalIteratorMethodTests(singletIterator);
        Assertions.assertEquals(list, list2);
        singletIterator.setAtom(testAtom1);
        Assertions.assertEquals(testAtom1, singletIterator.getAtom());
        list = generalIteratorMethodTests(singletIterator);
        Assertions.assertEquals(list, list1);
        Assertions.assertEquals(testAtom1, singletIterator.getAtom());
    }
    
    private AtomIteratorSinglet singletIterator;
    private IAtom testAtom1, testAtom2;
    private LinkedList list1, list2;

}
