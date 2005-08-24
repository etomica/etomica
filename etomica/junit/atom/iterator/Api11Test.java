package etomica.junit.atom.iterator;

import java.util.LinkedList;

import etomica.atom.Atom;
import etomica.atom.AtomPair;
import etomica.atom.AtomSet;
import etomica.atom.AtomTreeNodeGroup;
import etomica.atom.AtomsetArray;
import etomica.atom.SpeciesRoot;
import etomica.atom.iterator.Api11;
import etomica.junit.UnitTest;


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
public class Api11Test extends IteratorTest {

    /**
     * 
     */
    public Api11Test() {
        super();
//        UnitTest.VERBOSE = true;
    }
    
    public void setUp() {
        n0a = 5;
        nAtoms = 3;
        n1a = 10;
        n2a = 3;
        nTree = new int[] {5, 4, 3};
        SpeciesRoot root = UnitTest.makeStandardSpeciesTree(
                new int[] {n0a},
                nAtoms,
                new int[] {n1a},
                new int[] {n2a},
                nTree
        );
        rootNode = (AtomTreeNodeGroup)root.node;
        api11 = new Api11();
    }
    
    //species 0: 5 molecules, each a group of three atoms
    //species 1: 10 no-group single-atom molecules
    //species 2: 3 molecules, each with 5 subgroups of 4 groups of 3 atoms
    public void testIterator() {

        LinkedList list = null;
        AtomSet basis = null;
        Atom basis0 = null;
        Atom basis1 = null;
        AtomSet target = null;
        Atom target0 = null;
        Atom target1 = null;
        Atom target2 = null;
        Atom iterate0 = null;
        Atom iterate1 = null;
        
        assertEquals(api11.basisSize(), 2);
        
        //test initial iterator provides no iterates
        list = generalIteratorMethodTests(api11);
        assertEquals(list.size(), 0);
        
        //test targets different children of different basis atoms
        basis0 = rootNode.getDescendant(new int[] {0,0,0});
        target0 = rootNode.getDescendant(new int[] {0,0,0,1});
        basis1 = rootNode.getDescendant(new int[] {0,0,1});
        target1 = rootNode.getDescendant(new int[] {0,0,1,1});
        iterate0 = target0;
        iterate1 = target1;
        list = testOneIterate(basis0, basis1, target0, target1, iterate0, iterate1);

        //test targets different children of same basis atoms
        basis0 = rootNode.getDescendant(new int[] {0,0,0});
        target0 = rootNode.getDescendant(new int[] {0,0,0,0});
        basis1 = basis0;
        target1 = rootNode.getDescendant(new int[] {0,0,0,1});
        iterate0 = target0;
        iterate1 = target1;
        list = testOneIterate(basis0, basis1, target0, target1, iterate0, iterate1);

        //test targets different (grand)children of different basis atoms
        basis0 = rootNode.getDescendant(new int[] {0,2,0});
        target0 = rootNode.getDescendant(new int[] {0,2,0,1,1});
        basis1 = rootNode.getDescendant(new int[] {0,2,1});
        target1 = rootNode.getDescendant(new int[] {0,2,1,0});
        iterate0 = rootNode.getDescendant(new int[] {0,2,0,1});
        iterate1 = rootNode.getDescendant(new int[] {0,2,1,0});
        list = testOneIterate(basis0, basis1, target0, target1, iterate0, iterate1);

        //test both basis are different leaf atoms, and targets same as basis
        basis0 = rootNode.getDescendant(new int[] {0,0,0,1});
        target0 = basis0;
        basis1 = rootNode.getDescendant(new int[] {0,0,1,2});
        target1 = basis1;
        iterate0 = target0;
        iterate1 = target1;
        list = testOneIterate(basis0, basis1, target0, target1, iterate0, iterate1);

        //test one basis leaf atom, the other not, and targets same as basis
        basis0 = rootNode.getDescendant(new int[] {0,0,0,1});
        target0 = basis0;
        basis1 = rootNode.getDescendant(new int[] {0,0,1});
        target1 = basis1;
        testNoIterates(basis0, basis1, target0, target1);

        //test both basis are different leaf atoms, one target same, one target parent
        basis0 = rootNode.getDescendant(new int[] {0,0,0,1});
        target0 = basis0;
        basis1 = rootNode.getDescendant(new int[] {0,0,1,2});
        target1 = basis1.node.parentGroup();
        testNoIterates(basis0, basis1, target0, target1);

        //test target0, target1 both children of basis0, but not basis1 
        basis0 = rootNode.getDescendant(new int[] {0,0,0});
        target0 = rootNode.getDescendant(new int[] {0,0,0,1});
        basis1 = rootNode.getDescendant(new int[] {0,0,1});
        target1 = rootNode.getDescendant(new int[] {0,0,0,1});
        testNoIterates(basis0, basis1, target0, target1);

        //test neither target0, target1 children of basis0, basis1
        basis0 = rootNode.getDescendant(new int[] {0,0,0});
        target0 = rootNode.getDescendant(new int[] {0,1,0});
        basis1 = rootNode.getDescendant(new int[] {0,0,1});
        target1 = rootNode.getDescendant(new int[] {0,1,1});
        testNoIterates(basis0, basis1, target0, target1);
        
        //*************** test multiple targets
        
        //extra target doesn't specify an iterate
        basis0 = rootNode.getDescendant(new int[] {0,0,0});
        target0 = rootNode.getDescendant(new int[] {0,0,0,1});
        basis1 = rootNode.getDescendant(new int[] {0,0,1});
        target1 = rootNode.getDescendant(new int[] {0,0,1,1});
        target2 = rootNode.getDescendant(new int[] {0,1,0});
        testNoIterates(new AtomPair(basis0, basis1),
                new AtomsetArray(new Atom[] {target0, target1, target2}));

        //extra target specifies too many iterates
        basis0 = rootNode.getDescendant(new int[] {0,0,0});
        target0 = rootNode.getDescendant(new int[] {0,0,0,1});
        basis1 = rootNode.getDescendant(new int[] {0,0,1});
        target1 = rootNode.getDescendant(new int[] {0,0,1,1});
        target2 = rootNode.getDescendant(new int[] {0,0,0,2});//extra child of basis0
        testNoIterates(new AtomPair(basis0, basis1),
                new AtomsetArray(new Atom[] {target0, target1, target2}));

        //extra target specifies consistent iterate
        basis0 = rootNode.getDescendant(new int[] {0,2,0});
        target0 = rootNode.getDescendant(new int[] {0,2,0,1,0});
        basis1 = rootNode.getDescendant(new int[] {0,0,1});
        target1 = rootNode.getDescendant(new int[] {0,0,1,1});
        target2 = rootNode.getDescendant(new int[] {0,2,0,1,1});
        iterate0 = rootNode.getDescendant(new int[] {0,2,0,1});
        iterate1 = target1;
        testOneIterate(new AtomPair(basis0, basis1),
                new AtomsetArray(new Atom[] {target0, target1, target2}), 
                iterate0, iterate1);
        testOneIterate(new AtomPair(basis0, basis1),
                new AtomsetArray(new Atom[] {target1, target0, target2}), 
                iterate0, iterate1);
        testOneIterate(new AtomPair(basis0, basis1),
                new AtomsetArray(new Atom[] {target2, target0, target1}), 
                iterate0, iterate1);
        
        //test targets different conflicting children of same basis atoms
        basis0 = rootNode.getDescendant(new int[] {0,0,0});
        target0 = rootNode.getDescendant(new int[] {0,0,0,0});
        basis1 = basis0;
        target1 = rootNode.getDescendant(new int[] {0,0,0,1});
        target2 = rootNode.getDescendant(new int[] {0,0,0,2});
        iterate0 = target0;
        iterate1 = target1;
        testNoIterates(new AtomPair(basis0, basis1),
                new AtomsetArray(new Atom[] {target0, target1, target2}));

        //test targets different consistent children of same basis atoms
        basis0 = rootNode.getDescendant(new int[] {0,2,2});
        basis1 = basis0;
        target0 = rootNode.getDescendant(new int[] {0,2,2,0,0});
        target1 = rootNode.getDescendant(new int[] {0,2,2,1,0});
        target2 = rootNode.getDescendant(new int[] {0,2,2,1,1,1});
        iterate0 = rootNode.getDescendant(new int[] {0,2,2,0});
        iterate1 = rootNode.getDescendant(new int[] {0,2,2,1});
        testOneIterate(new AtomPair(basis0, basis1),
                new AtomsetArray(new Atom[] {target2, target0, target1}), 
                iterate0, iterate1);
        
        //test both basis are different leaf atoms, and targets same as basis
        basis0 = rootNode.getDescendant(new int[] {0,0,0,1});
        target0 = basis0;
        basis1 = rootNode.getDescendant(new int[] {0,0,1,2});
        target1 = basis1;
        target2 = target1;
        iterate0 = target0;
        iterate1 = target1;
        testOneIterate(new AtomPair(basis0, basis1),
                new AtomsetArray(new Atom[] {target0, target1, target2}), 
                iterate0, iterate1);

        //test for appropriate exceptions
        basis0 = rootNode.getDescendant(new int[] {0,0,0});
        basis1 = rootNode.getDescendant(new int[] {0,0,1});
        basis = new AtomsetArray(new Atom[] {basis0, basis1});
        target = basis;
//        AtomSet oneNull = new AtomsetArray(new Atom[] {basis0, null});
        testException(null, target);
        testException(AtomSet.NULL, target);
        testException(basis0, target);
        testException(basis, null);
        testException(basis, AtomSet.NULL);
        testException(basis, basis0);
//        testException(basis, oneNull);
//        testException(oneNull, target);
        

        
        //int[] {phase (0), species (0,1,2), molecule etc}
    }
    
    private LinkedList testOneIterate(Atom basis0, Atom basis1, 
            Atom target0, Atom target1, Atom iterate0, Atom iterate1) {
        LinkedList listA = testOneIterate(new AtomPair(basis0, basis1), 
                new AtomPair(target0, target1), iterate0, iterate1);
        LinkedList listB = testOneIterate(new AtomPair(basis0, basis1), 
                new AtomPair(target1, target0), iterate0, iterate1);
        assertEquals(listA, listB);
        return listA;
    }
    
    private LinkedList testOneIterate(AtomSet basis, AtomSet target, 
            Atom iterate0, Atom iterate1) {
        api11.setBasis(basis);
        assertTrue(api11.haveTarget(target));
        api11.setTarget(target);
        LinkedList list = generalIteratorMethodTests(api11);
        assertEquals(list.size(), 1);
        api11.reset();
        AtomSet pair = api11.next();
        if(basis.getAtom(0).equals(basis.getAtom(1))) {
            assertTrue(pair.getAtom(0).compareTo(pair.getAtom(1)) < 0);
        } else {
            assertTrue(pair.getAtom(0).node.isDescendedFrom(basis.getAtom(0)));
            assertTrue(pair.getAtom(1).node.isDescendedFrom(basis.getAtom(1)));
        }
        Lister testLister = new Lister();
        testLister.actionPerformed(new AtomPair(iterate0, iterate1));
        assertEquals(list, testLister.list);
        assertTrue(api11.haveTarget(target));//test again to ensure iteration didn't change anything
        return list;
    }
    
    private void testNoIterates(Atom basis0, Atom basis1, Atom target0, Atom target1) {
        testNoIterates(new AtomPair(basis0, basis1), new AtomPair(target0, target1));
        testNoIterates(new AtomPair(basis0, basis1), new AtomPair(target1, target0));
    }
    
    private void testNoIterates(AtomSet basis, AtomSet target) {
        api11.setBasis(basis);
        assertFalse(api11.haveTarget(target));
        api11.setTarget(target);
        LinkedList list = generalIteratorMethodTests(api11);
        assertEquals(list.size(), 0);
        assertFalse(api11.haveTarget(target));//test again to ensure iteration didn't change anything
    }
    
//    private void testException(Atom basis0, Atom basis1, AtomSet target) {
//        boolean exceptionThrown = false;
//        try {
//            api11.setBasis(new AtomsetArray(new Atom[] {basis0, basis1}));
//            api11.setTarget(target);
//        } catch(IllegalArgumentException ex) {
//            exceptionThrown = true;
//        }
//        assertTrue(exceptionThrown);
//    }
    
    private void testException(AtomSet basis, AtomSet target) {
        boolean exceptionThrown = false;
        try {
            api11.setBasis(basis);
            api11.setTarget(target);
        } catch(Exception ex) {
            exceptionThrown = true;
        }
        assertTrue(exceptionThrown);
    }
    
    private Api11 api11;
    private AtomTreeNodeGroup rootNode;
    int n0a, nAtoms, n1a, n2a;
    int[] nTree;

}
