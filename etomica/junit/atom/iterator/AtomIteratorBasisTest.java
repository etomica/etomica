package etomica.junit.atom.iterator;

import java.util.LinkedList;

import etomica.atom.Atom;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomList;
import etomica.atom.AtomSet;
import etomica.atom.AtomTreeNodeGroup;
import etomica.atom.AtomsetArray;
import etomica.atom.SpeciesRoot;
import etomica.atom.iterator.AtomIteratorBasis;
import etomica.junit.UnitTestUtil;


/**
 * Unit test for AtomIteratorSinglet class.
 *
 * @author David Kofke
 *
 */

/*
 * History
 * Created on May 30, 2005 by kofke
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
        n2a = 3;
        nTree = new int[] {5, 4, 3};
        SpeciesRoot root = UnitTestUtil.makeStandardSpeciesTree(
                new int[] {n0a},
                nAtoms,
                new int[] {n1a},
                new int[] {n2a},
                nTree
        );
        rootNode = (AtomTreeNodeGroup)root.node;
        
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
        Atom basis = null;
        Atom target = null;
        Atom iterate = null;
        AtomArrayList iterates = null;
        
        assertEquals(basisIterator.basisSize(), 1);
        
        //test initial iterator provides no iterates
        list = generalIteratorMethodTests(basisIterator);
        assertEquals(list.size(), 0);
        
        //test no-target iteration of children of a basis
        basis = rootNode.getDescendant(new int[] {0,0,0});
        target = null;
        iterates = (AtomArrayList)((AtomTreeNodeGroup)basis.node).childList.clone();
        list = testListIterates(basis, target, iterates);
        assertEquals(list.size(), nAtoms);

        //test no-target iteration of a leaf basis
        basis = rootNode.getDescendant(new int[] {0,0,0,1});
        target = null;
        iterate = basis;
        testOneIterate(basis, target, iterate);
        
        //test target is a child of the basis
        basis = rootNode.getDescendant(new int[] {0,0,0});
        target = rootNode.getDescendant(new int[] {0,0,0,1});
        iterate = (Atom)target;
        testOneIterate(basis, target, iterate);

        //test subsequent nulling of target
        basisIterator.setTarget(null);
        list = generalIteratorMethodTests(basisIterator);
        assertEquals(list.size(), nAtoms);
        testLister.clear();
        testLister.addEachToList(((AtomTreeNodeGroup)basis.node).childList.toArray());
        assertEquals(list, testLister.list);

        //test target is the basis, both not a leaf; should be same as target==null
        basis = rootNode.getDescendant(new int[] {0,0,0});
        target = basis;
        iterates = (AtomArrayList)((AtomTreeNodeGroup)basis.node).childList.clone();
        list = testListIterates(basis, target, iterates);
        assertEquals(list.size(), nAtoms);

        //test target is the basis, both a leaf
        basis = rootNode.getDescendant(new int[] {0,1,0});
        target = basis;
        iterate = basis;
        testOneIterate(basis, target, iterate);

        //test target is in hierarchy above basis, a leaf; should be same as target==null
        basis = rootNode.getDescendant(new int[] {0,0,0,0});
        target = rootNode.getDescendant(new int[] {0,0,0});
        iterate = basis;
        testOneIterate(basis, target, iterate);

        //test target is in hierarchy apart from basis; should return no iterates
        basis = rootNode.getDescendant(new int[] {0,0,0,0});
        target = rootNode.getDescendant(new int[] {0,1,0});
        testNoIterates(basis, target);

        //test target is derived from basis, but is not a child of it
        basis = rootNode.getDescendant(new int[] {0,2,0});
        target = rootNode.getDescendant(new int[] {0,2,0,1,0});
        iterate = rootNode.getDescendant(new int[] {0,2,0,1});
        testOneIterate(basis, target, iterate);
        
        //test specifying null target
        //also test specifying deeper basis
        basis = rootNode.getDescendant(new int[] {0,2,1,2});
        target = null;
        iterates = (AtomArrayList)((AtomTreeNodeGroup)basis.node).childList.clone();
        list = testListIterates(basis, target, iterates);
        
        //test null basis
        basis = null;
        target = rootNode.getDescendant(new int[] {0,2,0,1,0});
        testNoIterates(basis, target);
        
        //test null basis with null target
        basis = null;
        target = null;
        testNoIterates(basis, target);

        //int[] {phase (0), species (0,1,2), molecule etc}
    }
    
    private LinkedList testOneIterate(Atom basis, Atom target, Atom iterate) {
        basisIterator.setBasis(basis);
        assertTrue(basisIterator.haveTarget(target));
        basisIterator.setTarget(target);
        LinkedList list = generalIteratorMethodTests(basisIterator);
        Lister testLister = new Lister();
        testLister.actionPerformed(iterate);
        assertEquals(list, testLister.list);
        assertTrue(basisIterator.haveTarget(target));//test again to ensure iteration didn't change anything
        return list;
    }
    
    private LinkedList testListIterates(Atom basis, Atom target, AtomArrayList iterates) {
        basisIterator.setBasis(basis);
        assertTrue(basisIterator.haveTarget(target));
        basisIterator.setTarget(target);
        LinkedList list = generalIteratorMethodTests(basisIterator);
        Lister testLister = new Lister();
        testLister.addEachToList(iterates);
        assertEquals(list, testLister.list);
        assertTrue(basisIterator.haveTarget(target));//test again to ensure iteration didn't change anything
        return list;
    }
    
    private void testNoIterates(Atom basis, Atom target) {
        basisIterator.setBasis(basis);
        assertFalse(basisIterator.haveTarget(target));
        basisIterator.setTarget(target);
        LinkedList list = generalIteratorMethodTests(basisIterator);
        assertEquals(list.size(), 0);
        assertFalse(basisIterator.haveTarget(target));//test again to ensure iteration didn't change anything
    }
    
    private AtomIteratorBasis basisIterator;
    private AtomTreeNodeGroup rootNode;
    int n0a, nAtoms, n1a, n2a;
    int[] nTree;

}
