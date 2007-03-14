package etomica.junit.atom.iterator;

import java.util.LinkedList;

import etomica.atom.Atom;
import etomica.atom.AtomGroup;
import etomica.atom.AtomPair;
import etomica.atom.AtomType;
import etomica.atom.AtomTypeGroup;
import etomica.atom.SpeciesRoot;
import etomica.atom.iterator.ApiBuilder;
import etomica.atom.iterator.ApiIntergroup;
import etomica.atom.iterator.ApiIntragroup;
import etomica.atom.iterator.IteratorDirective;
import etomica.junit.UnitTestUtil;

/**
 * Tests the iterators made by the various static methods in ApiBuilder.
 * 
 * @author David Kofke
 *  
 */
public class ApiBuilderTest extends IteratorTestAbstract {

   public ApiBuilderTest() {
        super();
        UnitTestUtil.VERBOSE = false;
    }

   /**
    * Setup used by adjacentPairIterator and nonAdjacentPairIterator tests.
    *
    */
    public void setUpA() {
        n0a = 5;
        nAtoms = 10;
        n1a = 10;
        n2a = 3;
        nTree = new int[] { 5, 4, 3 };
        root = UnitTestUtil.makeStandardSpeciesTree(new int[] { n0a, 1 },
                nAtoms, new int[] { n1a, 2}, new int[] { n2a, 3 }, nTree);
    }
    
    /**
     * Sets up various choices of basis, target, direction and checks that
     * non-adjacent atoms are given in pairs with targeted iterate.
     */
    public void testNonAdjacentPairIterator() {
        setUpA();
        ApiIntragroup api = ApiBuilder.makeNonAdjacentPairIterator();

        setup1();
        nonAdjacentPairTests(api);
        
        setup2();
        api.setBasis(parent);
        
        //target matches basis, no direction, two iterates
        api.setDirection(null);
        api.setTarget(target);
        testApiIterates(api, target, upNon, dnNon);

        //target matches basis, up, one iterate
        api.setDirection(UP);
        testApiIterates(api, UP, target, upNon);
        
        //target matches basis, down, one iterate
        api.setDirection(DOWN);
        testApiIterates(api, DOWN, target, dnNon);
        
        //target doesn't match basis, no iterates
        api.setTarget(targetFirst);
        testNoIterates(api);
        
        setup3();
        nonAdjacentPairTests(api);

        setup4();
        nonAdjacentPairTests(api);

    }
    
    public void testIntergroupTypeIterator() {
        //make tree of two species
        //species 0 has 5 molecules, each with 5 atoms, 3 of one type, 2 of another
        //species 1 has 7 molecules, each with 11 atoms, 4 of one type, 1 of another, and 6 of another
        //iterator must loop over pairs formed from molecules of each species
        root = UnitTestUtil.makeMultitypeSpeciesTree(new int[] {5,7}, 
                new int[][] {{3,2},{4,1,6}});
        AtomTypeGroup rootType = (AtomTypeGroup)root.getType();
        AtomType[] types = new AtomType[2];
        AtomPair basisPair = new AtomPair();

        //test 3-atom type and 4-atom type, no target
        basisPair.atom0 = root.getDescendant(new int[] {0,0,2});
        basisPair.atom1 = root.getDescendant(new int[] {0,1,1});
        types[0] = rootType.getDescendant(new int[] {0,0,0,0});
        types[1] = rootType.getDescendant(new int[] {0,1,0,0});
        ApiIntergroup api = ApiBuilder.makeIntergroupTypeIterator(types);
        api.setBasis(basisPair);
        LinkedList list0 = generalIteratorMethodTests(api);
        assertEquals(list0.size(), 12);
        //test 3 and 4, one of the 3 given as target
        Atom target0 = root.getDescendant(new int[] {0,0,2,1});
        api.setTarget(target0);
        LinkedList list1 = generalIteratorMethodTests(api);
        assertEquals(list1.size(), 4);
        //test 3 and 4, one of the 4 given as target
        Atom target1 = root.getDescendant(new int[] {0,1,1,0});
        api.setTarget(target1);
        list1 = generalIteratorMethodTests(api);
        assertEquals(list1.size(), 3);
        //give target that isn't the specified type
        target0 = root.getDescendant(new int[] {0,0,2,4});
        api.setTarget(target0);
        testNoIterates(api);
        //again
        target1 = root.getDescendant(new int[] {0,1,1,10});
        api.setTarget(target1);
        testNoIterates(api);
        //no targets again
        api.setTarget(null);
        list1 = generalIteratorMethodTests(api);
        assertEquals(list0, list1);
        
        //same tests, but switch order of basis; nothing should give iterates
        //test 3-atom type and 4-atom type, no target
        basisPair.atom1 = root.getDescendant(new int[] {0,0,2});
        basisPair.atom0 = root.getDescendant(new int[] {0,1,1});
        types[0] = rootType.getDescendant(new int[] {0,0,0,0});
        types[1] = rootType.getDescendant(new int[] {0,1,0,0});
        api = ApiBuilder.makeIntergroupTypeIterator(types);
        api.setBasis(basisPair);
        testNoIterates(api);
        //test 3 and 4, one of the 3 given as target
        target0 = root.getDescendant(new int[] {0,0,2,1});
        api.setTarget(target0);
        testNoIterates(api);
        //test 3 and 4, one of the 4 given as target
        target1 = root.getDescendant(new int[] {0,1,1,0});
        api.setTarget(target1);
        testNoIterates(api);

        //same tests, but switch order of basis and switch order of types
        //test 3-atom type and 4-atom type, no target
        basisPair.atom1 = root.getDescendant(new int[] {0,0,2});
        basisPair.atom0 = root.getDescendant(new int[] {0,1,1});
        types[1] = rootType.getDescendant(new int[] {0,0,0,0});
        types[0] = rootType.getDescendant(new int[] {0,1,0,0});
        api = ApiBuilder.makeIntergroupTypeIterator(types);
        api.setBasis(basisPair);
        list0 = generalIteratorMethodTests(api);
        assertEquals(list0.size(), 12);
        //test 3 and 4, one of the 3 given as target
        target0 = root.getDescendant(new int[] {0,0,2,1});
        api.setTarget(target0);
        list1 = generalIteratorMethodTests(api);
        assertEquals(list1.size(), 4);
        //test 3 and 4, one of the 4 given as target
        target1 = root.getDescendant(new int[] {0,1,1,0});
        api.setTarget(target1);
        list1 = generalIteratorMethodTests(api);
        assertEquals(list1.size(), 3);

        //test 3-atom type and 1-atom type, no target
        basisPair.atom0 = root.getDescendant(new int[] {0,0,2});
        basisPair.atom1 = root.getDescendant(new int[] {0,1,1});
        types[0] = rootType.getDescendant(new int[] {0,0,0,0});
        types[1] = rootType.getDescendant(new int[] {0,1,0,1});
        api = ApiBuilder.makeIntergroupTypeIterator(types);
        api.setBasis(basisPair);
        list0 = generalIteratorMethodTests(api);
        assertEquals(list0.size(), 3);
        //test 3 and 1, one of the 3 given as target
        target0 = root.getDescendant(new int[] {0,0,2,1});
        api.setTarget(target0);
        list1 = generalIteratorMethodTests(api);
        assertEquals(list1.size(), 1);
        //test 3 and 1, the 1 given as target
        target1 = root.getDescendant(new int[] {0,1,1,4});
        api.setTarget(target1);
        list1 = generalIteratorMethodTests(api);
        assertEquals(list1.size(), 3);
        //give target that isn't the specified type
        target0 = root.getDescendant(new int[] {0,0,2,4});
        api.setTarget(target0);
        testNoIterates(api);
        //again
        target1 = root.getDescendant(new int[] {0,1,1,10});
        api.setTarget(target1);
        testNoIterates(api);
        //no targets again
        api.setTarget(null);
        list1 = generalIteratorMethodTests(api);
        assertEquals(list0, list1);

        basisPair.atom0 = root.getDescendant(new int[] {0,0,2});
        basisPair.atom1 = root.getDescendant(new int[] {0,1,1});
        types[0] = rootType.getDescendant(new int[] {0,0,0,0});
        types[1] = rootType.getDescendant(new int[] {0,1,0,0});
        api = ApiBuilder.makeIntergroupTypeIterator(types);
        api.setBasis(basisPair);
        list1 = generalIteratorMethodTests(api);
        assertEquals(list1.size(), 12);

        //incomplete
    }
    
    public void testIntragroupTypeIterator() {
        //AtomsetIteratorBasisDependent makeIntragroupTypeIterator(AtomType[] types)
        //incomplete
    }

    public void testAdjacentPairIterator() {
        setUpA();
        ApiIntragroup api = ApiBuilder.makeAdjacentPairIterator();
        
        setup1();
        adjacentPairTests(api);
        
        //************ test basis is leaf
        setup2();
        api.setBasis(parent);
        
        //target matches basis, no direction, two iterates
        api.setDirection(null);
        api.setTarget(target);
        testApiTwoIterates(api, new AtomPair(target, up), new AtomPair(dn, target));

        //target matches basis, up, one iterate
        api.setDirection(UP);
        testApiOneIterate(api, new AtomPair(target, up));
        
        //target matches basis, down, one iterate
        api.setDirection(DOWN);
        testApiOneIterate(api, new AtomPair(dn, target));
        
        //target doesn't match basis, no iterates
        api.setTarget(targetFirst);
        testNoIterates(api);

        //********* test target further descended from iterates
        setup3();

        adjacentPairTests(api);
        
        setup4();
        adjacentPairTests(api);

    }

    //******* adjacent/nonadjacent setup -- basis has only one child
    private void setup4() {
        parent = root.getDescendant(new int[] {1,0});//phase1, species0
        target = root.getDescendant(new int[] {1,0,0});//the only species0 molecule
        targetFirst = target;
        targetLast = target;
        up = dn = upFirst = dnLast = null;
        upNon = upFirstNon = dnNon = dnLastNon = new Atom[0]; 
        iterate = target;
        iterateFirst = target;
        iterateLast = target;
    }

    //************ adjacent/nonadjacent setup -- target is descended from but not direct child of basis
    private void setup3() {
        parent = root.getDescendant(new int[] {0,2,2});//phase0, species2, molecule2
        target = root.getDescendant(new int[] {0,2,2,1,0,1});
        targetFirst = root.getDescendant(new int[] {0,2,2,0,0,2});
        targetLast = root.getDescendant(new int[] {0,2,2,4,1});
        up = root.getDescendant(new int[] {0,2,2,2});
        upNon = new Atom[] {
                root.getDescendant(new int[] {0,2,2,3}),
                root.getDescendant(new int[] {0,2,2,4})};
        upFirst = root.getDescendant(new int[] {0,2,2,1});
        upFirstNon = new Atom[] {
                root.getDescendant(new int[] {0,2,2,2}),
                root.getDescendant(new int[] {0,2,2,3}),
                root.getDescendant(new int[] {0,2,2,4})};
        dn = root.getDescendant(new int[] {0,2,2,0});
        dnNon = new Atom[0];
        dnLast = root.getDescendant(new int[] {0,2,2,3});
        dnLastNon = new Atom[] {
                root.getDescendant(new int[] {0,2,2,2}),
                root.getDescendant(new int[] {0,2,2,1}),
                root.getDescendant(new int[] {0,2,2,0})};
        iterate = root.getDescendant(new int[] {0,2,2,1});
        iterateFirst = root.getDescendant(new int[] {0,2,2,0});
        iterateLast = root.getDescendant(new int[] {0,2,2,4});
    }


    //**********  adjacent/nonadjacent setup -- basis is a leaf atom
    private void setup2() {
        parent = root.getDescendant(new int[] {0,1,5});//leaf-atom basis
        target = root.getDescendant(new int[] {0,1,5});//atom5 
        targetFirst = root.getDescendant(new int[] {0,1,0});//atom0 
        targetLast = root.getDescendant(new int[] {0,1,9});//atom9
        up = root.getDescendant(new int[] {0,1,6});
        upNon = new Atom[] {
                root.getDescendant(new int[] {0,1,7}),
                root.getDescendant(new int[] {0,1,8}),
                root.getDescendant(new int[] {0,1,9})};
        upFirst = root.getDescendant(new int[] {0,1,1});
        upFirstNon = new Atom[] {
                root.getDescendant(new int[] {0,1,2}),
                root.getDescendant(new int[] {0,1,3}),
                root.getDescendant(new int[] {0,1,4}),
                root.getDescendant(new int[] {0,1,5}),
                root.getDescendant(new int[] {0,1,6}),
                root.getDescendant(new int[] {0,1,7}),
                root.getDescendant(new int[] {0,1,8}),
                root.getDescendant(new int[] {0,1,9})};
        dn = root.getDescendant(new int[] {0,1,4});
        dnNon = new Atom[] {
                root.getDescendant(new int[] {0,1,3}),
                root.getDescendant(new int[] {0,1,2}),
                root.getDescendant(new int[] {0,1,1}),
                root.getDescendant(new int[] {0,1,0})};
        dnLast = root.getDescendant(new int[] {0,1,8});
        dnLastNon = new Atom[] {
                root.getDescendant(new int[] {0,1,7}),
                root.getDescendant(new int[] {0,1,6}),
                root.getDescendant(new int[] {0,1,5}),
                root.getDescendant(new int[] {0,1,4}),
                root.getDescendant(new int[] {0,1,3}),
                root.getDescendant(new int[] {0,1,2}),
                root.getDescendant(new int[] {0,1,1}),
                root.getDescendant(new int[] {0,1,0})};

    }

    //******* adjacent/nonadjacent setup -- basis has child atoms, target is among them
    private void setup1() {
        parent = root.getDescendant(new int[] {0,0,2});
        target = root.getDescendant(new int[] {0,0,2,5});
        targetFirst = root.getDescendant(new int[] {0,0,2,0});
        targetLast = root.getDescendant(new int[] {0,0,2,9});
        up = root.getDescendant(new int[] {0,0,2,6});
        upNon = new Atom[] {
                        root.getDescendant(new int[] {0,0,2,7}),
                        root.getDescendant(new int[] {0,0,2,8}),
                        root.getDescendant(new int[] {0,0,2,9})};
        upFirst = root.getDescendant(new int[] {0,0,2,1});
        upFirstNon = new Atom[] {
                        root.getDescendant(new int[] {0,0,2,2}),
                        root.getDescendant(new int[] {0,0,2,3}),
                        root.getDescendant(new int[] {0,0,2,4}),
                        root.getDescendant(new int[] {0,0,2,5}),
                        root.getDescendant(new int[] {0,0,2,6}),
                        root.getDescendant(new int[] {0,0,2,7}),
                        root.getDescendant(new int[] {0,0,2,8}),
                        root.getDescendant(new int[] {0,0,2,9})};
        dn = root.getDescendant(new int[] {0,0,2,4});
        dnNon = new Atom[] {
                        root.getDescendant(new int[] {0,0,2,3}),
                        root.getDescendant(new int[] {0,0,2,2}),
                        root.getDescendant(new int[] {0,0,2,1}),
                        root.getDescendant(new int[] {0,0,2,0})};
        dnLast = root.getDescendant(new int[] {0,0,2,8});
        dnLastNon = new Atom[] {
                        root.getDescendant(new int[] {0,0,2,7}),
                        root.getDescendant(new int[] {0,0,2,6}),
                        root.getDescendant(new int[] {0,0,2,5}),
                        root.getDescendant(new int[] {0,0,2,4}),
                        root.getDescendant(new int[] {0,0,2,3}),
                        root.getDescendant(new int[] {0,0,2,2}),
                        root.getDescendant(new int[] {0,0,2,1}),
                        root.getDescendant(new int[] {0,0,2,0})};
        iterate = target;
        iterateFirst = targetFirst;
        iterateLast = targetLast;
    }

    /**
     * Test AdjacentPairIterator.
     * @param api the iterator
     */
    private void adjacentPairTests(ApiIntragroup api) {
        
        api.setBasis(parent);

        //target, no direction, two iterates
        api.setDirection(null);
        api.setTarget(target);
        testApiTwoIterates(api, new AtomPair(iterate, up), new AtomPair(dn, iterate));
        
        //first in list, no direction, one iterate
        api.setDirection(null);
        api.setTarget(targetFirst);
        testApiOneIterate(api, new AtomPair(iterateFirst, upFirst));
        
        //first in list, down, no iterates
        api.setDirection(DOWN);
        api.setTarget(targetFirst);
        testNoIterates(api);

        //last in list, down, one iterate
        api.setDirection(DOWN);
        api.setTarget(targetLast);
        testApiOneIterate(api, new AtomPair(dnLast, iterateLast));

        //last in list, no direction, one iterate
        api.setDirection(null);
        api.setTarget(targetLast);
        testApiOneIterate(api, new AtomPair(dnLast, iterateLast));
        
        //last in list, up, no iterates
        api.setDirection(UP);
        api.setTarget(targetLast);
        testNoIterates(api);

        //first in list, up, one iterate
        api.setDirection(UP);
        api.setTarget(targetFirst);
        testApiOneIterate(api, new AtomPair(iterateFirst, upFirst));


        //target, up, one iterate
        api.setDirection(UP);
        api.setTarget(target);
        testApiOneIterate(api, new AtomPair(iterate, up));

        //target, down, one iterate
        api.setDirection(DOWN);
        api.setTarget(target);
        testApiOneIterate(api, new AtomPair(dn, iterate));

        //no target, n-1 iterates
        api.setTarget(null);
        LinkedList list0 = generalIteratorMethodTests(api);
        assertEquals(list0.size(), ((AtomGroup)parent).getChildList().size()-1);
        
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
     */
    private void nonAdjacentPairTests(ApiIntragroup api) {
        
        api.setBasis(parent);

        //target, no direction, two iterates
        api.setDirection(null);
        api.setTarget(target);
        testApiIterates(api, iterate, upNon, dnNon);
        
        //first in list, no direction, one iterate
        api.setDirection(null);
        api.setTarget(targetFirst);
        testApiIterates(api, UP, iterateFirst, upFirstNon);
        
        //first in list, down, no iterates
        api.setDirection(DOWN);
        api.setTarget(targetFirst);
        testNoIterates(api);

        //last in list, down, one iterate
        api.setDirection(DOWN);
        api.setTarget(targetLast);
        testApiIterates(api, DOWN, iterateLast, dnLastNon);

        //last in list, no direction, one iterate
        api.setDirection(null);
        api.setTarget(targetLast);
        testApiIterates(api, DOWN, iterateLast, dnLastNon);
        
        //last in list, up, no iterates
        api.setDirection(UP);
        api.setTarget(targetLast);
        testNoIterates(api);

        //first in list, up, one iterate
        api.setDirection(UP);
        api.setTarget(targetFirst);
        testApiIterates(api, UP, iterateFirst, upFirstNon);


        //target, up, one iterate
        api.setDirection(UP);
        api.setTarget(target);
        testApiIterates(api, UP, iterate, upNon);

        //target, down, one iterate
        api.setDirection(DOWN);
        api.setTarget(target);
        testApiIterates(api, DOWN, iterate, dnNon);

        //no target, (2(n-2) + (n-2)*(n-3))/2 iterates
        //(first and last have n-2, the other n-2 have n-3)
        api.setTarget(null);
        LinkedList list0 = generalIteratorMethodTests(api);
        int n = ((AtomGroup)parent).getChildList().size();
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

    private SpeciesRoot root;
    int n0a, nAtoms, n1a, n2a, n3a;
    int[] nTree;
    private static final IteratorDirective.Direction UP = IteratorDirective.Direction.UP;
    private static final IteratorDirective.Direction DOWN = IteratorDirective.Direction.DOWN;
    private Atom parent;
    private Atom target;
    private Atom targetFirst;
    private Atom targetLast;
    private Atom up;
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
