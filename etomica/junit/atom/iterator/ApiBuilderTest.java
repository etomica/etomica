package etomica.junit.atom.iterator;

import java.util.LinkedList;

import etomica.atom.AtomArrayList;
import etomica.atom.AtomGroup;
import etomica.atom.AtomPair;
import etomica.atom.AtomType;
import etomica.atom.AtomTypeGroup;
import etomica.atom.IAtom;
import etomica.atom.IAtomGroup;
import etomica.atom.SpeciesAgent;
import etomica.atom.AtomManager;
import etomica.atom.iterator.ApiBuilder;
import etomica.atom.iterator.ApiIntergroup;
import etomica.atom.iterator.ApiIntragroup;
import etomica.atom.iterator.IteratorDirective;
import etomica.junit.UnitTestUtil;
import etomica.simulation.Simulation;

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
        sim = UnitTestUtil.makeStandardSpeciesTree(new int[] { n0a, 1 },
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
        sim = UnitTestUtil.makeMultitypeSpeciesTree(new int[] {5,7}, 
                new int[][] {{3,2},{4,1,6}});
        AtomType[] types = new AtomType[2];
        AtomPair basisPair = new AtomPair();

        AtomManager atomManager = sim.getPhases()[0].getSpeciesMaster();
        SpeciesAgent agent0 = (SpeciesAgent)atomManager.getAgentList().get(0);
        SpeciesAgent agent1 = (SpeciesAgent)atomManager.getAgentList().get(1);
        AtomTypeGroup agentType0 = sim.getSpeciesManager().getSpeciesAgentTypes()[0];
        AtomTypeGroup agentType1 = sim.getSpeciesManager().getSpeciesAgentTypes()[1];
        
        //test 3-atom type and 4-atom type, no target
        basisPair.atom0 = agent0.getDescendant(new int[] {2});
        basisPair.atom1 = agent1.getDescendant(new int[] {1});
        types[0] = agentType0.getDescendant(new int[] {0,0});
        types[1] = agentType1.getDescendant(new int[] {0,0});
        ApiIntergroup api = ApiBuilder.makeIntergroupTypeIterator(types);
        api.setBasis(basisPair);
        LinkedList list0 = generalIteratorMethodTests(api);
        assertEquals(list0.size(), 12);
        //test 3 and 4, one of the 3 given as target
        IAtom target0 = agent0.getDescendant(new int[] {2,1});
        api.setTarget(target0);
        LinkedList list1 = generalIteratorMethodTests(api);
        assertEquals(list1.size(), 4);
        //test 3 and 4, one of the 4 given as target
        IAtom target1 = agent1.getDescendant(new int[] {1,0});
        api.setTarget(target1);
        list1 = generalIteratorMethodTests(api);
        assertEquals(list1.size(), 3);
        //give target that isn't the specified type
        target0 = agent0.getDescendant(new int[] {2,4});
        api.setTarget(target0);
        testNoIterates(api);
        //again
        target1 = agent1.getDescendant(new int[] {1,10});
        api.setTarget(target1);
        testNoIterates(api);
        //no targets again
        api.setTarget(null);
        list1 = generalIteratorMethodTests(api);
        assertEquals(list0, list1);
        
        //same tests, but switch order of basis; nothing should give iterates
        //test 3-atom type and 4-atom type, no target
        basisPair.atom1 = agent0.getDescendant(new int[] {2});
        basisPair.atom0 = agent1.getDescendant(new int[] {1});
        types[0] = agentType0.getDescendant(new int[] {0,0});
        types[1] = agentType1.getDescendant(new int[] {0,0});
        api = ApiBuilder.makeIntergroupTypeIterator(types);
        api.setBasis(basisPair);
        testNoIterates(api);
        //test 3 and 4, one of the 3 given as target
        target0 = agent0.getDescendant(new int[] {2,1});
        api.setTarget(target0);
        testNoIterates(api);
        //test 3 and 4, one of the 4 given as target
        target1 = agent1.getDescendant(new int[] {1,0});
        api.setTarget(target1);
        testNoIterates(api);

        //same tests, but switch order of basis and switch order of types
        //test 3-atom type and 4-atom type, no target
        basisPair.atom1 = agent0.getDescendant(new int[] {2});
        basisPair.atom0 = agent1.getDescendant(new int[] {1});
        types[1] = agentType0.getDescendant(new int[] {0,0});
        types[0] = agentType1.getDescendant(new int[] {0,0});
        api = ApiBuilder.makeIntergroupTypeIterator(types);
        api.setBasis(basisPair);
        list0 = generalIteratorMethodTests(api);
        assertEquals(list0.size(), 12);
        //test 3 and 4, one of the 3 given as target
        target0 = agent0.getDescendant(new int[] {2,1});
        api.setTarget(target0);
        list1 = generalIteratorMethodTests(api);
        assertEquals(list1.size(), 4);
        //test 3 and 4, one of the 4 given as target
        target1 = agent1.getDescendant(new int[] {1,0});
        api.setTarget(target1);
        list1 = generalIteratorMethodTests(api);
        assertEquals(list1.size(), 3);

        //test 3-atom type and 1-atom type, no target
        basisPair.atom0 = agent0.getDescendant(new int[] {2});
        basisPair.atom1 = agent1.getDescendant(new int[] {1});
        types[0] = agentType0.getDescendant(new int[] {0,0});
        types[1] = agentType1.getDescendant(new int[] {0,1});
        api = ApiBuilder.makeIntergroupTypeIterator(types);
        api.setBasis(basisPair);
        list0 = generalIteratorMethodTests(api);
        assertEquals(list0.size(), 3);
        //test 3 and 1, one of the 3 given as target
        target0 = agent0.getDescendant(new int[] {2,1});
        api.setTarget(target0);
        list1 = generalIteratorMethodTests(api);
        assertEquals(list1.size(), 1);
        //test 3 and 1, the 1 given as target
        target1 = agent1.getDescendant(new int[] {1,4});
        api.setTarget(target1);
        list1 = generalIteratorMethodTests(api);
        assertEquals(list1.size(), 3);
        //give target that isn't the specified type
        target0 = agent0.getDescendant(new int[] {2,4});
        api.setTarget(target0);
        testNoIterates(api);
        //again
        target1 = agent1.getDescendant(new int[] {1,10});
        api.setTarget(target1);
        testNoIterates(api);
        //no targets again
        api.setTarget(null);
        list1 = generalIteratorMethodTests(api);
        assertEquals(list0, list1);

        basisPair.atom0 = agent0.getDescendant(new int[] {2});
        basisPair.atom1 = agent1.getDescendant(new int[] {1});
        types[0] = agentType0.getDescendant(new int[] {0,0});
        types[1] = agentType1.getDescendant(new int[] {0,0});
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
        AtomManager atomManager = sim.getPhases()[1].getSpeciesMaster();
        parent = atomManager.getAgentList().get(0);//phase1, species0
        target = ((IAtomGroup)parent).getChildList().get(0);//the only species0 molecule
        targetFirst = target;
        targetLast = target;
        up = dn = upFirst = dnLast = null;
        upNon = upFirstNon = dnNon = dnLastNon = new IAtom[0]; 
        iterate = target;
        iterateFirst = target;
        iterateLast = target;
    }

    //************ adjacent/nonadjacent setup -- target is descended from but not direct child of basis
    private void setup3() {
        AtomManager atomManager = sim.getPhases()[0].getSpeciesMaster();
        parent = ((IAtomGroup)atomManager.getAgentList().get(2)).getChildList().get(2); //Descendant(new int[] {2,2});//phase0, species2, molecule2
        AtomArrayList childList = ((IAtomGroup)parent).getChildList(); 
        target = ((AtomGroup)parent).getDescendant(new int[] {1,0,1});
        targetFirst = ((AtomGroup)parent).getDescendant(new int[] {0,0,2});
        targetLast = ((AtomGroup)parent).getDescendant(new int[] {4,1});
        up = ((AtomGroup)parent).getDescendant(new int[] {2});
        upNon = new IAtom[] {childList.get(3), childList.get(4)};
        upFirst = childList.get(1);
        upFirstNon = new IAtom[] {childList.get(2), childList.get(3),
                childList.get(4)};
        dn = childList.get(0);
        dnNon = new IAtom[0];
        dnLast = childList.get(3);
        dnLastNon = new IAtom[] {childList.get(2),childList.get(1),
                childList.get(0)};
        iterate = childList.get(1);
        iterateFirst = childList.get(0);
        iterateLast = childList.get(4);
    }


    //**********  adjacent/nonadjacent setup -- basis is a leaf atom
    private void setup2() {
        AtomManager atomManager = sim.getPhases()[0].getSpeciesMaster();
        AtomArrayList moleculeList = ((IAtomGroup)atomManager.getAgentList().get(1)).getChildList();
        parent = moleculeList.get(5);//leaf-atom basis
        target = parent;//atom5 
        targetFirst = moleculeList.get(0);//atom0 
        targetLast = moleculeList.get(9);//atom9
        up = moleculeList.get(6);
        upNon = new IAtom[] {moleculeList.get(7),moleculeList.get(8),
                moleculeList.get(9)};
        upFirst = moleculeList.get(1);
        upFirstNon = new IAtom[] {moleculeList.get(2),moleculeList.get(3),
                moleculeList.get(4),moleculeList.get(5),moleculeList.get(6),
                moleculeList.get(7),moleculeList.get(8),moleculeList.get(9)};
        dn = moleculeList.get(4);
        dnNon = new IAtom[] {moleculeList.get(3),moleculeList.get(2),
                moleculeList.get(1),moleculeList.get(0)};
        dnLast = moleculeList.get(8);
        dnLastNon = new IAtom[] {moleculeList.get(7),moleculeList.get(6),
                moleculeList.get(5),moleculeList.get(4),moleculeList.get(3),
                moleculeList.get(2),moleculeList.get(1),moleculeList.get(0)};
    }

    //******* adjacent/nonadjacent setup -- basis has child atoms, target is among them
    private void setup1() {
        AtomManager atomManager = sim.getPhases()[0].getSpeciesMaster();
        parent = ((IAtomGroup)atomManager.getAgentList().get(0)).getChildList().get(2);
        AtomArrayList childList = ((IAtomGroup)parent).getChildList();
        target = childList.get(5);
        targetFirst = childList.get(0);
        targetLast = childList.get(9);
        up = childList.get(6);
        upNon = new IAtom[] {childList.get(7),childList.get(8),childList.get(9)};
        upFirst = childList.get(1);
        upFirstNon = new IAtom[] {childList.get(2),childList.get(3),
                childList.get(4),childList.get(5),childList.get(6),
                childList.get(7),childList.get(8),childList.get(9)};
        dn = childList.get(4);
        dnNon = new IAtom[] {childList.get(3),childList.get(2),childList.get(1),
                childList.get(0)};
        dnLast = childList.get(8);
        dnLastNon = new IAtom[] {childList.get(7),childList.get(6),
                childList.get(5),childList.get(4),childList.get(3),
                childList.get(2),childList.get(1),childList.get(0)};
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
        assertEquals(list0.size(), ((IAtomGroup)parent).getChildList().size()-1);
        
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
        int n = ((IAtomGroup)parent).getChildList().size();
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

    private Simulation sim;
    int n0a, nAtoms, n1a, n2a, n3a;
    int[] nTree;
    private static final IteratorDirective.Direction UP = IteratorDirective.Direction.UP;
    private static final IteratorDirective.Direction DOWN = IteratorDirective.Direction.DOWN;
    private IAtom parent;
    private IAtom target;
    private IAtom targetFirst;
    private IAtom targetLast;
    private IAtom up;
    private IAtom[] upNon;
    private IAtom upFirst;
    private IAtom[] upFirstNon;
    private IAtom dn;
    private IAtom[] dnNon;
    private IAtom dnLast;
    private IAtom[] dnLastNon;
    private IAtom iterate;
    private IAtom iterateFirst;
    private IAtom iterateLast;
}
