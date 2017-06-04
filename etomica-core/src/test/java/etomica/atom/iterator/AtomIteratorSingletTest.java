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
        testAtom1 = new Atom(Space3D.getInstance());
        testAtom2 = new Atom(Space3D.getInstance());
        list1 = makeTestList(new IAtomList[] {new AtomSetSinglet(testAtom1)});
        list2 = makeTestList(new IAtomList[] {new AtomSetSinglet(testAtom2)});
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
