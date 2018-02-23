/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom.iterator;

import etomica.UnitTestUtil;
import etomica.atom.MoleculesetAction;
import etomica.molecule.*;
import etomica.molecule.iterator.MoleculeIterator;
import etomica.molecule.iterator.MoleculesetIterator;
import etomica.potential.IteratorDirective;
import junit.framework.TestCase;

import java.util.LinkedList;

/**
 * Provides library methods to test the basic functioning of an iterator.
 * Most JUnit iterator tests extend this class.
 * 
 * @author Ken Benjamin and David Kofke
 *  
 */
public abstract class MoleculeIteratorTestAbstract extends TestCase {


    /**
     * Clears the list of each ListerMolecule in the array.
     */
    public void clearLists(ListerMolecule[] ListerMolecule) {
        for (int i = 0; i < ListerMolecule.length; i++) {
            ListerMolecule[i].list.clear();
        }
    }

    /**
     * Prints all of the lists.
     */
    public void printLists(ListerMolecule[] ListerMolecule) {
        if (UnitTestUtil.VERBOSE) {
            for (int i = 0; i < ListerMolecule.length; i++) {
                System.out.println(ListerMolecule[i]);
            }
            System.out.println();
        }
    }
    
    /**
     * Returns a list from the given atoms, which can be checked against
     * the list returned by generalIteratorMethodTests.
     */
    public LinkedList makeTestList(IMoleculeList[] atoms) {
        ListerMolecule ListerMolecule = new ListerMolecule();
        for(int i=0; i<atoms.length; i++) {
            ListerMolecule.actionPerformed(atoms[i]);
        }
        return ListerMolecule.list;
    }
    
    protected void print(String string) {
        if(UnitTestUtil.VERBOSE) System.out.println(string);
    }

    public static void allAtoms(MoleculesetIterator mpi, MoleculesetAction action) {
        mpi.reset();
        for (IMoleculeList molecules = mpi.next(); molecules != null; molecules = mpi.next()) {
            action.actionPerformed(molecules);
        }
    }

    /**
     * Tests all the elementary methods of an atomset iterator. Methods tested are
     * allAtoms, hasNext/next, peek, size, contains, nBody, unset, reset.  Suitable
     * for use by any AtomsetIterator.  Tests nextAtom if iterator is instance
     * of AtomIterator
     * 
     * @param iterator
     *            an iterator in a condition to be tested; reset and iteration
     *            will be performed repeatedly on the iterator
     * @return a list of ListerMolecule instances for the atoms given by the iterator
     */
    protected java.util.LinkedList<String> generalIteratorMethodTests(
            MoleculesetIterator iterator) {
        ListerMolecule[] listerMolecule = ListerMolecule.listerArray(2);

        //******* test of next
        iterator.reset();
        for (IMoleculeList atomSet = iterator.next(); atomSet != null;
             atomSet = iterator.next()) {
            listerMolecule[0].actionPerformed(atomSet);
        }
        //******* test of size
        print("Testing size -- iterator.size, ListerMolecule.list.size: "+iterator.size()+" "+listerMolecule[0].list.size());
        assertEquals(iterator.size(), listerMolecule[0].list.size());

        IMoleculeList[] atoms = new IMoleculeList[listerMolecule[0].list.size()];
        
        //******* test of next
        iterator.reset();
        int j=0;
        for (IMoleculeList nextAtom = iterator.next(); nextAtom != null;
               nextAtom = iterator.next()) {
            atoms[j] = new MoleculesetArray(nextAtom);
            listerMolecule[1].actionPerformed(nextAtom);
            j++;
        }

        //******* test of nextAtom
        if(iterator instanceof MoleculeIterator) {
            iterator.reset();
            int i = 0;
            for (IMolecule next = ((MoleculeIterator)iterator).nextMolecule(); i<atoms.length;
                 next = ((MoleculeIterator)iterator).nextMolecule()) {
                assertEquals(next, atoms[i++].get(0));
            }
            assertNull(((MoleculeIterator)iterator).nextMolecule());
        }
        
        //******* test of unset
        iterator.reset();
        iterator.unset();
        assertNull(iterator.next());
        //******* test that next calls cannot cause next to return not-null
        for (int i = 0; i < 5; i++) {
            assertNull(iterator.next());
        }
        print("Just tested unset method");

        //******* test of nBody
        iterator.reset();
        for (IMoleculeList atomSet = iterator.next(); atomSet != null;
             atomSet = iterator.next()) {
            assertEquals(atomSet.size(), iterator.nBody());
        }

        print("Just tested nBody method");
        
//        if(iterator.nBody() == 2) {
//            iterator.reset();
//            while(iterator.hasNext()) {
//                AtomSet pair = iterator.next();
//                System.out.println(pair.getAtom(0).node.index() + " " + pair.getAtom(1).node.index());
//                assertTrue(pair.getAtom(0).node.index() < pair.getAtom(1).node.index());
//            }
//        }

        return listerMolecule[0].list;
    }

    /**
     * Tests that iterator gives no iterates
     */
    protected void testNoIterates(MoleculesetIterator iterator) {
        LinkedList list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(), 0);
    }

    /**
     * Tests that iterates given by iterator match a specified list.
     * @param iterator the conditioned iterator
     * @param partners array of expected iterates
     * @return the ListerMolecule list of iterates
     */
    protected LinkedList testIterates(MoleculeIterator iterator, IMolecule[] iterates) {
        LinkedList list = generalIteratorMethodTests(iterator);
        ListerMolecule test = new ListerMolecule();
        for(int i=0; i<iterates.length; i++) {
            test.actionPerformed(new MoleculeSetSinglet(iterates[i]));
        }
        assertEquals(list, test.list);
        return list;
    }

    /**
     * Tests that iterates given by iterate match a specified list, for
     * iteration uplist or dnlist (not both)
     * @param iterator the conditioned iterate
     * @param iterate the atom0 expected in all pair iterates
     * @param partners array of atom1 expected in the pair iterates
     * @return the ListerMolecule list of iterates
     */
    protected LinkedList testApiIterates(MoleculesetIterator iterator, IteratorDirective.Direction direction,
            IMolecule iterate, IMolecule[] partners) {
        LinkedList list = generalIteratorMethodTests(iterator);
        ListerMolecule test = new ListerMolecule();
        for(int i=0; i<partners.length; i++) {
            if(direction == IteratorDirective.Direction.UP) {
                test.actionPerformed(new MoleculePair(iterate, partners[i]));
            } else {
                test.actionPerformed(new MoleculePair(partners[i], iterate));
            }
        }
        assertEquals(list, test.list);
        return list;
    }
    
    /**
     * Same as testApiIterates, but with atom1 the same in all pair iterates, while
     * atom0 varies.
     */
    protected LinkedList testApiIteratesSwap(MoleculesetIterator iterator, IMolecule iterate, IMolecule[] partners) {
        LinkedList list = generalIteratorMethodTests(iterator);
        ListerMolecule test = new ListerMolecule();
        for(int i=0; i<partners.length; i++) {
            test.actionPerformed(new MoleculePair(partners[i],iterate));
        }
        assertEquals(list, test.list);
        return list;
    }


    /**
     * Tests that iterator gives two particular iterates
     */
    protected LinkedList testApiTwoIterates(MoleculesetIterator iterator, MoleculePair pair0, MoleculePair pair1) {
        if((pair0.mol1 == null || pair0.mol0 == null) && (pair1.mol0 == null || pair1.mol1 == null)) {
            testNoIterates(iterator);
            return new LinkedList();
        }
        LinkedList list = generalIteratorMethodTests(iterator);
        ListerMolecule test = new ListerMolecule();
        test.actionPerformed(pair0);
        test.actionPerformed(pair1);
        assertEquals(list, test.list);
        return list;
    }

    /**
     * Tests that iterator gives a single particular iterate.
     */
    protected LinkedList testApiOneIterate(MoleculesetIterator iterator, MoleculePair pair) {
        if(pair.mol0 == null || pair.mol1 == null) {
            testNoIterates(iterator);
            return new LinkedList();
        }
        LinkedList list = generalIteratorMethodTests(iterator);
        ListerMolecule test = new ListerMolecule();
        test.actionPerformed(pair);
        assertEquals(list, test.list);
        return list;
    }

    /**
     * Tests that iterator gives a specific number of iterates.
     * @param api the iterator
     * @param n the number of iterates it should give
     * @return the ListerMolecule list of iterates
     */
    protected LinkedList countTest(MoleculesetIterator api, int n) {
        LinkedList list = generalIteratorMethodTests(api);
        assertEquals(list.size(), n);
        return list;
    }

    /**
     * Tests that iterates given by iterate match a specified list, for
     * iteration uplist and then dnlist
     * @param iterator the conditioned iterate
     * @param iterate the atom0 expected in all pair iterates
     * @param up array of atom1 expected in the uplist pair iterates
     * @param dn array of atom1 expected in the dnlist pair iterates
     * @return the ListerMolecule list of iterates
     */
    protected LinkedList testApiIterates(MoleculesetIterator iterator, IMolecule iterate, IMolecule[] up, IMolecule[] dn) {
        LinkedList list = generalIteratorMethodTests(iterator);
        ListerMolecule test = new ListerMolecule();
        for(int i=0; i<up.length; i++) {
            test.actionPerformed(new MoleculePair(iterate, up[i]));
        }
        for(int i=0; i<dn.length; i++) {
            test.actionPerformed(new MoleculePair(dn[i], iterate));
        }
        assertEquals(list, test.list);
        return list;
    }

}
