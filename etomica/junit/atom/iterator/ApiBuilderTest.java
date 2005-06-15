package etomica.junit.atom.iterator;

import java.util.LinkedList;

import etomica.Atom;
import etomica.AtomPair;
import etomica.AtomPairIterator;
import etomica.AtomSet;
import etomica.AtomTreeNodeGroup;
import etomica.AtomType;
import etomica.AtomsetIterator;
import etomica.IteratorDirective;
import etomica.SpeciesRoot;
import etomica.atom.iterator.ApiBuilder;
import etomica.atom.iterator.ApiIntergroup;
import etomica.atom.iterator.ApiIntragroup;
import etomica.atom.iterator.AtomsetIteratorBasisDependent;
import etomica.junit.UnitTest;

/**
 * Tests the iterators made by the various methods in ApiBuilder.
 * 
 * @author David Kofke
 *  
 */

/*
 * History Created on Jun 5, 2005 by kofke
 */
public class ApiBuilderTest extends IteratorTest {

    /**
     *  
     */
    public ApiBuilderTest() {
        super();
        UnitTest.VERBOSE = false;
    }

    public void setUp() {
        n0a = 5;
        nAtoms = 10;
        n1a = 10;
        n2a = 3;
        nTree = new int[] { 5, 4, 3 };
        SpeciesRoot root = UnitTest.makeStandardSpeciesTree(new int[] { n0a },
                nAtoms, new int[] { n1a }, new int[] { n2a }, nTree);
        rootNode = (AtomTreeNodeGroup) root.node;
    }
    
    public void testNonAdjacentPairIterator() {
        ApiIntragroup api = ApiBuilder.makeNonAdjacentPairIterator();
        //incomplete
    }
    
    public void testIntergroupTypeIterator() {
        //ApiIntergroup api makeIntergroupTypeIterator(AtomType[] types)
        //incomplete
    }
    
    public void testIntragroupTypeIterator() {
        //AtomsetIteratorBasisDependent makeIntragroupTypeIterator(AtomType[] types)
        //incomplete
    }

    public void testAdjacentPairIterator() {
        
        ApiIntragroup api = ApiBuilder.makeAdjacentPairIterator();
        
        //********** test targets are iterates
        Atom parent = rootNode.getDescendant(new int[] {0,0,2});//phase0, species0, molecule2
        Atom target = rootNode.getDescendant(new int[] {0,0,2,5});//atom5 of molecule
        Atom targetFirst = rootNode.getDescendant(new int[] {0,0,2,0});//atom0 of molecule
        Atom targetLast = rootNode.getDescendant(new int[] {0,0,2,9});//atom9 of molecule
        Atom up = rootNode.getDescendant(new int[] {0,0,2,6});
        Atom upFirst = rootNode.getDescendant(new int[] {0,0,2,1});
        Atom dn = rootNode.getDescendant(new int[] {0,0,2,4});
        Atom dnLast = rootNode.getDescendant(new int[] {0,0,2,8});
        Atom iterate = target;
        Atom iterateFirst = targetFirst;
        Atom iterateLast = targetLast;
        
        adjacentPairTests(api, parent, target, targetFirst, targetLast, 
                iterate, iterateFirst, iterateLast,
                up, upFirst, dn, dnLast);
        
        //************ test basis is leaf
        parent = rootNode.getDescendant(new int[] {0,1,5});//leaf-atom basis
        target = rootNode.getDescendant(new int[] {0,1,5});//atom5 
        targetFirst = rootNode.getDescendant(new int[] {0,1,0});//atom0 
        targetLast = rootNode.getDescendant(new int[] {0,1,9});//atom9
        up = rootNode.getDescendant(new int[] {0,1,6});
        upFirst = rootNode.getDescendant(new int[] {0,1,1});
        dn = rootNode.getDescendant(new int[] {0,1,4});
        dnLast = rootNode.getDescendant(new int[] {0,1,8});
        api.setBasis(parent);
        
        //target matches basis, no direction, two iterates
        api.setDirection(null);
        api.setTarget(target);
        testTwoIterates(api, new AtomPair(target, up), new AtomPair(target, dn));

        //target matches basis, up, one iterate
        api.setDirection(UP);
        testOneIterate(api, new AtomPair(target, up));
        
        //target matches basis, down, one iterate
        api.setDirection(DOWN);
        testOneIterate(api, new AtomPair(target, dn));
        
        //target doesn't match basis, no iterates
        api.setTarget(targetFirst);
        testNoIterates(api);

        //********* test target further descended from iterates
        parent = rootNode.getDescendant(new int[] {0,2,2});//phase0, species2, molecule2
        target = rootNode.getDescendant(new int[] {0,2,2,1,0,1});
        targetFirst = rootNode.getDescendant(new int[] {0,2,2,0,0,2});
        targetLast = rootNode.getDescendant(new int[] {0,2,2,4,1});
        up = rootNode.getDescendant(new int[] {0,2,2,2});
        upFirst = rootNode.getDescendant(new int[] {0,2,2,1});
        dn = rootNode.getDescendant(new int[] {0,2,2,0});
        dnLast = rootNode.getDescendant(new int[] {0,2,2,3});
        iterate = rootNode.getDescendant(new int[] {0,2,2,1});
        iterateFirst = rootNode.getDescendant(new int[] {0,2,2,0});
        iterateLast = rootNode.getDescendant(new int[] {0,2,2,4});

        adjacentPairTests(api, parent, target, targetFirst, targetLast,
                iterate, iterateFirst, iterateLast,
               up, upFirst, dn, dnLast);
    }

    private void adjacentPairTests(ApiIntragroup api, Atom parent, 
            Atom target, Atom targetFirst, Atom targetLast,
            Atom iterate, Atom iterateFirst, Atom iterateLast,
            Atom up, Atom upFirst, 
            Atom dn, Atom dnLast) {
        
        api.setBasis(parent);

        //target, no direction, two iterates
        api.setDirection(null);
        api.setTarget(target);
        testTwoIterates(api, new AtomPair(iterate, up), new AtomPair(iterate, dn));
        
        //first in list, no direction, one iterate
        api.setDirection(null);
        api.setTarget(targetFirst);
        testOneIterate(api, new AtomPair(iterateFirst, upFirst));
        
        //first in list, down, no iterates
        api.setDirection(DOWN);
        api.setTarget(targetFirst);
        testNoIterates(api);

        //last in list, down, one iterate
        api.setDirection(DOWN);
        api.setTarget(targetLast);
        testOneIterate(api, new AtomPair(iterateLast, dnLast));

        //last in list, no direction, one iterate
        api.setDirection(null);
        api.setTarget(targetLast);
        testOneIterate(api, new AtomPair(iterateLast, dnLast));
        
        //last in list, up, no iterates
        api.setDirection(UP);
        api.setTarget(targetLast);
        testNoIterates(api);

        //first in list, up, one iterate
        api.setDirection(UP);
        api.setTarget(targetFirst);
        testOneIterate(api, new AtomPair(iterateFirst, upFirst));


        //target, up, one iterate
        api.setDirection(UP);
        api.setTarget(target);
        testOneIterate(api, new AtomPair(iterate, up));

        //target, down, one iterate
        api.setDirection(DOWN);
        api.setTarget(target);
        testOneIterate(api, new AtomPair(iterate, dn));

        //no target, n-1 iterates
        api.setTarget(AtomSet.NULL);
        LinkedList list0 = generalIteratorMethodTests(api);
        assertEquals(list0.size(), ((AtomTreeNodeGroup)parent.node).childList.size()-1);
        
        //if no target, direction doesn't matter
        api.setDirection(null);
        LinkedList list1 = generalIteratorMethodTests(api);
        assertEquals(list0, list1);

        //no basis, no iterates
        api.setBasis(null);
        testNoIterates(api);
        api.setTarget(target);
        testNoIterates(api);
    }

    private LinkedList testTwoIterates(AtomPairIterator iterator, AtomPair pair0, AtomPair pair1) {
        LinkedList list = generalIteratorMethodTests(iterator);
        Lister test = new Lister();
        test.actionPerformed(pair0);
        test.actionPerformed(pair1);
        assertEquals(list, test.list);
        return list;
    }

    private LinkedList testOneIterate(AtomPairIterator iterator, AtomPair pair) {
        LinkedList list = generalIteratorMethodTests(iterator);
        Lister test = new Lister();
        test.actionPerformed(pair);
        assertEquals(list, test.list);
        return list;
    }
    
    private void testNoIterates(AtomsetIterator iterator) {
        LinkedList list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(), 0);
    }

    
    private AtomTreeNodeGroup rootNode;
    int n0a, nAtoms, n1a, n2a;
    int[] nTree;
    private final IteratorDirective.Direction UP = IteratorDirective.UP;
    private final IteratorDirective.Direction DOWN = IteratorDirective.DOWN;


}