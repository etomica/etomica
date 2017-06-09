/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom.iterator;

import java.util.LinkedList;

import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.atom.IMolecule;
import etomica.atom.IMoleculeList;
import etomica.simulation.Simulation;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomSetSinglet;
import etomica.atom.MoleculeSetSinglet;
import etomica.UnitTestUtil;


/**
 * Unit test for AtomIteratorSinglet class.
 *
 * @author David Kofke
 *
 */
public class AtomIteratorBasisTest extends IteratorTestAbstract {

    /**
     * 
     */
    public AtomIteratorBasisTest() {
        super();
    }
    
    public void setUp() {
        n0a = 5;
        nAtoms = 3;
        n1a = 10;
        sim = UnitTestUtil.makeStandardSpeciesTree(
                new int[] {n0a},
                nAtoms,
                new int[] {n1a}
        );
        
//        root = UnitTest.makeStandardSpeciesTree(
//                new int[] {n0a, n0b},
//                nAtoms,
//                new int[] {n1a, n1b},
//                new int[] {n2a, n2b},
//                nTree
//        );
        basisIterator = new AtomIteratorBasis();
    }
    
    //species 0: 5 molecules, each a group of three atoms
    //species 1: 10 no-group single-atom molecules
    //species 2: 3 molecules, each with 5 subgroups of 4 groups of 3 atoms
    public void testIterator() {

        Lister testLister = new Lister();
        LinkedList list = null;
        IMolecule basis = null;
        IAtom target = null;
        IAtom iterate = null;
        AtomArrayList iterates = null;
        
        assertEquals(basisIterator.basisSize(), 1);
        
        //test initial iterator provides no iterates
        list = generalIteratorMethodTests(basisIterator);
        assertEquals(list.size(), 0);
        
        //test no-target iteration of children of a basis
        Box box = sim.getBox(0);
        IMoleculeList moleculeList0 = box.getMoleculeList(sim.getSpecies(0));
        basis = moleculeList0.getMolecule(0);
        target = null;
        iterates = new AtomArrayList();
        iterates.addAll(basis.getChildList());
        list = testListIterates(basis, target, iterates);
        assertEquals(list.size(), nAtoms);

//        //test no-target iteration of a leaf basis
//        basis = ((IMolecule)moleculeList0.getAtom(0)).getChildList().getAtom(1);
//        target = null;
//        iterate = basis;
//        testOneIterate(basis, target, iterate);
        
        //test target is a child of the basis
        basis = moleculeList0.getMolecule(0);
        target = moleculeList0.getMolecule(0).getChildList().getAtom(1);
        iterate = target;
        testOneIterate(basis, target, iterate);

        //test subsequent nulling of target
        basisIterator.setTarget(null);
        list = generalIteratorMethodTests(basisIterator);
        assertEquals(list.size(), nAtoms);
        testLister.clear();
        testLister.addEachToList(basis.getChildList());
        assertEquals(list, testLister.list);

        //test specifying null target
        basis = moleculeList0.getMolecule(0);
        target = null;
        iterates.clear();
        iterates.addAll(basis.getChildList());
        list = testListIterates(basis, target, iterates);
        
        //test null basis
        basis = null;
        target = moleculeList0.getMolecule(0).getChildList().getAtom(1);
        testNoIterates(basis, target);
        
        //test null basis with null target
        basis = null;
        target = null;
        testNoIterates(basis, target);

        //int[] {box (0), species (0,1,2), molecule etc}
    }
    
    private LinkedList testOneIterate(IMolecule basis, IAtom target, IAtom iterate) {
        basisIterator.setBasis(new MoleculeSetSinglet(basis));
        assertTrue(basisIterator.haveTarget(target));
        basisIterator.setTarget(target);
        LinkedList list = generalIteratorMethodTests(basisIterator);
        Lister testLister = new Lister();
        testLister.actionPerformed(new AtomSetSinglet(iterate));
        assertEquals(list, testLister.list);
        assertTrue(basisIterator.haveTarget(target));//test again to ensure iteration didn't change anything
        return list;
    }
    
    private LinkedList testListIterates(IMolecule basis, IAtom target, AtomArrayList iterates) {
        basisIterator.setBasis(new MoleculeSetSinglet(basis));
        assertTrue(basisIterator.haveTarget(target));
        basisIterator.setTarget(target);
        LinkedList list = generalIteratorMethodTests(basisIterator);
        Lister testLister = new Lister();
        testLister.addEachToList(iterates);
        assertEquals(list, testLister.list);
        assertTrue(basisIterator.haveTarget(target));//test again to ensure iteration didn't change anything
        return list;
    }
    
    private void testNoIterates(IMolecule basis, IAtom target) {
        basisIterator.setBasis(basis == null ? null : new MoleculeSetSinglet(basis));
        assertFalse(basisIterator.haveTarget(target));
        basisIterator.setTarget(target);
        LinkedList list = generalIteratorMethodTests(basisIterator);
        assertEquals(list.size(), 0);
        assertFalse(basisIterator.haveTarget(target));//test again to ensure iteration didn't change anything
    }
    
    private AtomIteratorBasis basisIterator;
    private Simulation sim;
    int n0a, nAtoms, n1a;

}
