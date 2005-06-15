package etomica.junit.atom.iterator;

import java.util.LinkedList;

import etomica.Atom;
import etomica.AtomPair;
import etomica.AtomPairIterator;
import etomica.AtomSet;
import etomica.AtomTreeNodeGroup;
import etomica.AtomsetIterator;
import etomica.IteratorDirective;
import etomica.SpeciesRoot;
import etomica.atom.iterator.ApiBuilder;
import etomica.atom.iterator.ApiIntragroup;
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
        SpeciesRoot root = UnitTest.makeStandardSpeciesTree(new int[] { n0a, 1 },
                nAtoms, new int[] { n1a, 2}, new int[] { n2a, 3 }, nTree);
        rootNode = (AtomTreeNodeGroup) root.node;
    }
    
    /**
     * Sets up various choices of basis, target, direction and checks that
     * non-adjacent atoms are given in pairs with targeted iterate.
     */
    public void testNonAdjacentPairIterator() {
        ApiIntragroup api = ApiBuilder.makeNonAdjacentPairIterator();

        setup1();
        nonAdjacentPairTests(api, parent, target, targetFirst, targetLast, 
                iterate, iterateFirst, iterateLast,
                upNon, upFirstNon, dnNon, dnLastNon);
        
        setup2();
        api.setBasis(parent);
        
        //target matches basis, no direction, two iterates
        api.setDirection(null);
        api.setTarget(target);
        testIterates(api, target, upNon, dnNon);

        //target matches basis, up, one iterate
        api.setDirection(UP);
        testIterates(api, target, upNon);
        
        //target matches basis, down, one iterate
        api.setDirection(DOWN);
        testIterates(api, target, dnNon);
        
        //target doesn't match basis, no iterates
        api.setTarget(targetFirst);
        testNoIterates(api);
        
        setup3();
        nonAdjacentPairTests(api, parent, target, targetFirst, targetLast, 
                iterate, iterateFirst, iterateLast,
                upNon, upFirstNon, dnNon, dnLastNon);

        setup4();
        nonAdjacentPairTests(api, parent, target, targetFirst, targetLast, 
                iterate, iterateFirst, iterateLast,
                upNon, upFirstNon, dnNon, dnLastNon);

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
        
        setup1();
        adjacentPairTests(api, parent, target, targetFirst, targetLast, 
                iterate, iterateFirst, iterateLast,
                up2, upFirst, dn, dnLast);
        
        //************ test basis is leaf
        setup2();
        api.setBasis(parent);
        
        //target matches basis, no direction, two iterates
        api.setDirection(null);
        api.setTarget(target);
        testTwoIterates(api, new AtomPair(target, up2), new AtomPair(target, dn));

        //target matches basis, up, one iterate
        api.setDirection(UP);
        testOneIterate(api, new AtomPair(target, up2));
        
        //target matches basis, down, one iterate
        api.setDirection(DOWN);
        testOneIterate(api, new AtomPair(target, dn));
        
        //target doesn't match basis, no iterates
        api.setTarget(targetFirst);
        testNoIterates(api);

        //********* test target further descended from iterates
        setup3();

        adjacentPairTests(api, parent, target, targetFirst, targetLast,
                iterate, iterateFirst, iterateLast,
               up2, upFirst, dn, dnLast);
        
        setup4();
        adjacentPairTests(api, parent, target, targetFirst, targetLast,
                iterate, iterateFirst, iterateLast,
               up2, upFirst, dn, dnLast);

    }

    //******* basis has only one child
    private void setup4() {
        parent = rootNode.getDescendant(new int[] {1,0});//phase1, species0
        target = rootNode.getDescendant(new int[] {1,0,0});
        targetFirst = target;
        targetLast = target;
        up2 = dn = upFirst = dnLast = null;
        upNon = upFirstNon = dnNon = dnLastNon = new Atom[0]; 
        iterate = target;
        iterateFirst = target;
        iterateLast = target;
    }

    //************ target is descended from but not direct child of basis
    private void setup3() {
        parent = rootNode.getDescendant(new int[] {0,2,2});//phase0, species2, molecule2
        target = rootNode.getDescendant(new int[] {0,2,2,1,0,1});
        targetFirst = rootNode.getDescendant(new int[] {0,2,2,0,0,2});
        targetLast = rootNode.getDescendant(new int[] {0,2,2,4,1});
        up2 = rootNode.getDescendant(new int[] {0,2,2,2});
        upNon = new Atom[] {
                rootNode.getDescendant(new int[] {0,2,2,3}),
                rootNode.getDescendant(new int[] {0,2,2,4})};
        upFirst = rootNode.getDescendant(new int[] {0,2,2,1});
        upFirstNon = new Atom[] {
                rootNode.getDescendant(new int[] {0,2,2,2}),
                rootNode.getDescendant(new int[] {0,2,2,3}),
                rootNode.getDescendant(new int[] {0,2,2,4})};
        dn = rootNode.getDescendant(new int[] {0,2,2,0});
        dnNon = new Atom[0];
        dnLast = rootNode.getDescendant(new int[] {0,2,2,3});
        dnLastNon = new Atom[] {
                rootNode.getDescendant(new int[] {0,2,2,2}),
                rootNode.getDescendant(new int[] {0,2,2,1}),
                rootNode.getDescendant(new int[] {0,2,2,0})};
        iterate = rootNode.getDescendant(new int[] {0,2,2,1});
        iterateFirst = rootNode.getDescendant(new int[] {0,2,2,0});
        iterateLast = rootNode.getDescendant(new int[] {0,2,2,4});
    }


    //**********  basis is a leaf atom
    private void setup2() {
        parent = rootNode.getDescendant(new int[] {0,1,5});//leaf-atom basis
        target = rootNode.getDescendant(new int[] {0,1,5});//atom5 
        targetFirst = rootNode.getDescendant(new int[] {0,1,0});//atom0 
        targetLast = rootNode.getDescendant(new int[] {0,1,9});//atom9
        up2 = rootNode.getDescendant(new int[] {0,1,6});
        upNon = new Atom[] {
                rootNode.getDescendant(new int[] {0,1,7}),
                rootNode.getDescendant(new int[] {0,1,8}),
                rootNode.getDescendant(new int[] {0,1,9})};
        upFirst = rootNode.getDescendant(new int[] {0,1,1});
        upFirstNon = new Atom[] {
                rootNode.getDescendant(new int[] {0,1,2}),
                rootNode.getDescendant(new int[] {0,1,3}),
                rootNode.getDescendant(new int[] {0,1,4}),
                rootNode.getDescendant(new int[] {0,1,5}),
                rootNode.getDescendant(new int[] {0,1,6}),
                rootNode.getDescendant(new int[] {0,1,7}),
                rootNode.getDescendant(new int[] {0,1,8}),
                rootNode.getDescendant(new int[] {0,1,9})};
        dn = rootNode.getDescendant(new int[] {0,1,4});
        dnNon = new Atom[] {
                rootNode.getDescendant(new int[] {0,1,3}),
                rootNode.getDescendant(new int[] {0,1,2}),
                rootNode.getDescendant(new int[] {0,1,1}),
                rootNode.getDescendant(new int[] {0,1,0})};
        dnLast = rootNode.getDescendant(new int[] {0,1,8});
        dnLastNon = new Atom[] {
                rootNode.getDescendant(new int[] {0,1,7}),
                rootNode.getDescendant(new int[] {0,1,6}),
                rootNode.getDescendant(new int[] {0,1,5}),
                rootNode.getDescendant(new int[] {0,1,4}),
                rootNode.getDescendant(new int[] {0,1,3}),
                rootNode.getDescendant(new int[] {0,1,2}),
                rootNode.getDescendant(new int[] {0,1,1}),
                rootNode.getDescendant(new int[] {0,1,0})};

    }

    //******* basis has child atoms, target is among them
    private void setup1() {
        parent = rootNode.getDescendant(new int[] {0,0,2});
        target = rootNode.getDescendant(new int[] {0,0,2,5});
        targetFirst = rootNode.getDescendant(new int[] {0,0,2,0});
        targetLast = rootNode.getDescendant(new int[] {0,0,2,9});
        up2 = rootNode.getDescendant(new int[] {0,0,2,6});
        upNon = new Atom[] {
                        rootNode.getDescendant(new int[] {0,0,2,7}),
                        rootNode.getDescendant(new int[] {0,0,2,8}),
                        rootNode.getDescendant(new int[] {0,0,2,9})};
        upFirst = rootNode.getDescendant(new int[] {0,0,2,1});
        upFirstNon = new Atom[] {
                        rootNode.getDescendant(new int[] {0,0,2,2}),
                        rootNode.getDescendant(new int[] {0,0,2,3}),
                        rootNode.getDescendant(new int[] {0,0,2,4}),
                        rootNode.getDescendant(new int[] {0,0,2,5}),
                        rootNode.getDescendant(new int[] {0,0,2,6}),
                        rootNode.getDescendant(new int[] {0,0,2,7}),
                        rootNode.getDescendant(new int[] {0,0,2,8}),
                        rootNode.getDescendant(new int[] {0,0,2,9})};
        dn = rootNode.getDescendant(new int[] {0,0,2,4});
        dnNon = new Atom[] {
                        rootNode.getDescendant(new int[] {0,0,2,3}),
                        rootNode.getDescendant(new int[] {0,0,2,2}),
                        rootNode.getDescendant(new int[] {0,0,2,1}),
                        rootNode.getDescendant(new int[] {0,0,2,0})};
        dnLast = rootNode.getDescendant(new int[] {0,0,2,8});
        dnLastNon = new Atom[] {
                        rootNode.getDescendant(new int[] {0,0,2,7}),
                        rootNode.getDescendant(new int[] {0,0,2,6}),
                        rootNode.getDescendant(new int[] {0,0,2,5}),
                        rootNode.getDescendant(new int[] {0,0,2,4}),
                        rootNode.getDescendant(new int[] {0,0,2,3}),
                        rootNode.getDescendant(new int[] {0,0,2,2}),
                        rootNode.getDescendant(new int[] {0,0,2,1}),
                        rootNode.getDescendant(new int[] {0,0,2,0})};
        iterate = target;
        iterateFirst = targetFirst;
        iterateLast = targetLast;
    }

    /**
     * Test AdjacentPairIterator.
     * @param api the iterator
     * @param parent the basis atom
     * @param target a target in middle of basis list
     * @param targetFirst a target at beginning of basis list
     * @param targetLast a target at end of basis list
     * @param iterate the child-of-basis iterate implied by the target (same as target if target is in childlist of basis)
     * @param iterateFirst the child-of-basis iterate implied by targetFirst
     * @param iterateLast the child-of-basis iterate implied by targetLast
     * @param up the adjacent atom uplist of iterate
     * @param upFirst the adjacent atom uplist of iterateFirst
     * @param dn the adjacent atom dnlist of iterate
     * @param dnLast the adjacent atom dnlist of iterateLast
     */
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
    
    /**
     * Test NonAdjacentPairIterator.
     * @param api the iterator
     * @param parent the basis atom
     * @param target a target in middle of basis list
     * @param targetFirst a target at beginning of basis list
     * @param targetLast a target at end of basis list
     * @param iterate the child-of-basis iterate implied by the target (same as target if target is in childlist of basis)
     * @param iterateFirst the child-of-basis iterate implied by targetFirst
     * @param iterateLast the child-of-basis iterate implied by targetLast
     * @param up the nonadjacent atoms uplist of iterate
     * @param upFirst the nonadjacent atoms uplist of iterateFirst
     * @param dn the nonadjacent atoms dnlist of iterate
     * @param dnLast the nonadjacent atoms dnlist of iterateLast
     */

    private void nonAdjacentPairTests(ApiIntragroup api, Atom parent, 
            Atom target, Atom targetFirst, Atom targetLast,
            Atom iterate, Atom iterateFirst, Atom iterateLast,
            Atom[] up, Atom[] upFirst, 
            Atom[] dn, Atom[] dnLast) {
        
        api.setBasis(parent);

        //target, no direction, two iterates
        api.setDirection(null);
        api.setTarget(target);
        testIterates(api, iterate, up, dn);
        
        //first in list, no direction, one iterate
        api.setDirection(null);
        api.setTarget(targetFirst);
        testIterates(api, iterateFirst, upFirst);
        
        //first in list, down, no iterates
        api.setDirection(DOWN);
        api.setTarget(targetFirst);
        testNoIterates(api);

        //last in list, down, one iterate
        api.setDirection(DOWN);
        api.setTarget(targetLast);
        testIterates(api, iterateLast, dnLast);

        //last in list, no direction, one iterate
        api.setDirection(null);
        api.setTarget(targetLast);
        testIterates(api, iterateLast, dnLast);
        
        //last in list, up, no iterates
        api.setDirection(UP);
        api.setTarget(targetLast);
        testNoIterates(api);

        //first in list, up, one iterate
        api.setDirection(UP);
        api.setTarget(targetFirst);
        testIterates(api, iterateFirst, upFirst);


        //target, up, one iterate
        api.setDirection(UP);
        api.setTarget(target);
        testIterates(api, iterate, up);

        //target, down, one iterate
        api.setDirection(DOWN);
        api.setTarget(target);
        testIterates(api, iterate, dn);

        //no target, (2(n-2) + (n-2)*(n-3))/2 iterates
        //(first and last have n-2, the other n-2 have n-3)
        api.setTarget(AtomSet.NULL);
        LinkedList list0 = generalIteratorMethodTests(api);
        int n = ((AtomTreeNodeGroup)parent.node).childList.size();
        assertEquals(list0.size(), (2*(n-2) + (n-2)*(n-3))/2);
        
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

    /**
     * Tests that iterates given by iterate match a specified list, for
     * iteration uplist and then dnlist
     * @param iterator the conditioned iterate
     * @param iterate the atom0 expected in all pair iterates
     * @param up array of atom1 expected in the uplist pair iterates
     * @param dn array of atom1 expected in the dnlist pair iterates
     * @return the Lister list of iterates
     */
    private LinkedList testIterates(AtomPairIterator iterator, Atom iterate,
            Atom[] up, Atom[] dn) {
        LinkedList list = generalIteratorMethodTests(iterator);
        Lister test = new Lister();
        for(int i=0; i<up.length; i++) {
            test.actionPerformed(new AtomPair(iterate, up[i]));
        }
        for(int i=0; i<dn.length; i++) {
            test.actionPerformed(new AtomPair(iterate, dn[i]));
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
     * @return the Lister list of iterates
     */
    private LinkedList testIterates(AtomPairIterator iterator, Atom iterate,
            Atom[] partners) {
        LinkedList list = generalIteratorMethodTests(iterator);
        Lister test = new Lister();
        for(int i=0; i<partners.length; i++) {
            test.actionPerformed(new AtomPair(iterate, partners[i]));
        }
        assertEquals(list, test.list);
        return list;
    }

    /**
     * Tests that iterator gives two particular iterates
     */
    private LinkedList testTwoIterates(AtomPairIterator iterator, AtomPair pair0, AtomPair pair1) {
        if(pair0.atom1 == null && pair1.atom1 == null) {
            testNoIterates(iterator);
            return new LinkedList();
        }
        LinkedList list = generalIteratorMethodTests(iterator);
        Lister test = new Lister();
        test.actionPerformed(pair0);
        test.actionPerformed(pair1);
        assertEquals(list, test.list);
        return list;
    }

    /**
     * Tests that iterator gives a single particular iterate.
     */
    private LinkedList testOneIterate(AtomPairIterator iterator, AtomPair pair) {
        if(pair.atom1 == null) {
            testNoIterates(iterator);
            return new LinkedList();
        }
        LinkedList list = generalIteratorMethodTests(iterator);
        Lister test = new Lister();
        test.actionPerformed(pair);
        assertEquals(list, test.list);
        return list;
    }

    /**
     * Tests that iterator gives no iterates
     */
    private void testNoIterates(AtomsetIterator iterator) {
        LinkedList list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(), 0);
    }

    
    private AtomTreeNodeGroup rootNode;
    int n0a, nAtoms, n1a, n2a, n3a;
    int[] nTree;
    private final IteratorDirective.Direction UP = IteratorDirective.UP;
    private final IteratorDirective.Direction DOWN = IteratorDirective.DOWN;
    private Atom parent;
    private Atom target;
    private Atom targetFirst;
    private Atom targetLast;
    private Atom up2;
    private Atom[] upNon;
    private Atom upFirst;
    private Atom[] upFirstNon;
    private Atom dn;
    private Atom[] dnNon;
    private Atom dnLast;
    private Atom[] dnLastNon;
    private Atom iterate;
    private Atom iterateFirst;
    private Atom iterateLast;


}