package etomica.junit.atom.iterator;

import java.util.LinkedList;

import etomica.atom.AtomArrayList;
import etomica.atom.AtomSet;
import etomica.atom.IAtom;
import etomica.atom.SpeciesAgent;
import etomica.atom.iterator.AtomIteratorTreePhase;
import etomica.junit.UnitTestUtil;
import etomica.phase.Phase;
import etomica.simulation.Simulation;
import etomica.species.Species;

/**
 * Unit test for AtomIteratorTree
 */
public class AtomIteratorTreePhaseTest extends IteratorTestAbstract {

    protected void setUp() {
        n0a = 3;
        nAtoms = 3;
        n1a = 10;
        n2a = 4;
        nTree = new int[] { 5, 4, 3 };
        Simulation sim = UnitTestUtil.makeStandardSpeciesTree(new int[] { n0a },
                nAtoms, new int[] { n1a }, new int[] { n2a }, nTree);
        phase = sim.getPhases()[0];
        Species[] species = sim.getSpeciesManager().getSpecies();
        speciesAgent0 = phase.getAgent(species[0]);
        speciesAgent1 = phase.getAgent(species[1]);
        speciesAgent2 = phase.getAgent(species[2]);

        treeIterator = new AtomIteratorTreePhase();
        treeIterator.setPhase(phase);
    }
    
    //species 0: 5 molecules, each a group of three atoms
    //species 1: 10 no-group single-atom molecules
    //species 2: 3 molecules, each with 5 subgroups of 4 groups of 3 atoms
    public void testIterator() {
        LinkedList list = null;
        int count = 0;

        //test iteration over all nodes from speciesMaster
        count = 3 + n0a*(1 + nAtoms) + n1a*(1) 
                + n2a*(1 + nTree[0]*(1 + nTree[1]*(1 + nTree[2])));
        treeIterator.setDoAllNodes(true);
        list = testIterateCount(count);

        //test iteration over different depths, starting at root
        treeIterator.setDoAllNodes(true);

        list = testArrayIterates(1, new IAtom[] {speciesAgent0, speciesAgent1, speciesAgent2});
        list = testIterateCount(2, 3+n0a+n1a+n2a);
        list = testIterateCount(3, 3+n0a*(1+nAtoms)+n1a+n2a*(1+nTree[0]));
        list = testIterateCount(4, 3+n0a*(1+nAtoms)+n1a+n2a*(1+nTree[0]*(1+nTree[1])));
        list = testIterateCount(5, 3+n0a*(1+nAtoms)+n1a+n2a*(1+nTree[0]*(1+nTree[1]*(1+nTree[2]))));
        list = testIterateCount(100, 3+n0a*(1+nAtoms)+n1a+n2a*(1+nTree[0]*(1+nTree[1]*(1+nTree[2]))));

        treeIterator.setDoAllNodes(false);
        list = testArrayIterates(1, new IAtom[] {speciesAgent0, speciesAgent1, speciesAgent2});
        list = testIterateCount(2, n0a+n1a+n2a);
        list = testIterateCount(3, n0a*nAtoms+n1a+n2a*nTree[0]);
        list = testIterateCount(4, n0a*nAtoms+n1a+n2a*nTree[0]*(nTree[1]));
        list = testIterateCount(5, n0a*nAtoms+n1a+n2a*nTree[0]*nTree[1]*nTree[2]);
        list = testIterateCount(100, n0a*nAtoms+n1a+n2a*nTree[0]*nTree[1]*nTree[2]);

        //test iteration of leaf atoms
        treeIterator.setAsLeafIterator();
        count = n0a*nAtoms + n1a + n2a*nTree[0]*nTree[1]*nTree[2];
        list = testIterateCount(count);
        list = testListIterates(phase.getSpeciesMaster().getLeafList());
        
        //test re-specifying iteration in different orders
        int depth = 3;
        boolean doAllNodes = true;
        count = 3+n0a*(1+nAtoms)+n1a+n2a*(1+nTree[0]);
        treeIterator = new AtomIteratorTreePhase(phase, depth, doAllNodes);
        list = generalIteratorMethodTests(treeIterator);
        assertEquals(list.size(), count);
        treeIterator = new AtomIteratorTreePhase();
        treeIterator.setDoAllNodes(doAllNodes);
        treeIterator.setPhase(phase);
        treeIterator.setIterationDepth(depth);
        LinkedList list1 = generalIteratorMethodTests(treeIterator);
        assertEquals(list, list1);

        treeIterator.setIterationDepth(depth);
        treeIterator.setDoAllNodes(doAllNodes);
        treeIterator.setPhase(phase);
        list1 = generalIteratorMethodTests(treeIterator);
        assertEquals(list, list1);
        
        treeIterator.setPhase(phase);
        treeIterator.setIterationDepth(depth);
        treeIterator.setDoAllNodes(doAllNodes);
        list1 = generalIteratorMethodTests(treeIterator);
        assertEquals(list, list1);
    }
    
    private LinkedList testIterateCount(int depth, int count) {
        treeIterator.setIterationDepth(depth);
        LinkedList list0 = testIterateCount(count);
        treeIterator.setIterationDepth(1);
        treeIterator.setIterationDepth(depth);
        LinkedList list1 = generalIteratorMethodTests(treeIterator);
        assertEquals(list0, list1);
        return list0;
    }
    private LinkedList testIterateCount(int count) {
        LinkedList list = generalIteratorMethodTests(treeIterator);
//        System.out.println(list.size()+ " "+count);
        assertEquals(list.size(), count);
        return list;
    }
    
    private LinkedList testArrayIterates(int depth, IAtom[] iterates) {
        treeIterator.setIterationDepth(depth);
        AtomArrayList arrayList = new AtomArrayList();
        for (int i=0; i<iterates.length; i++) {
            arrayList.add(iterates[i]);
        }
        return testListIterates(arrayList);
    }
    
    private LinkedList testListIterates(AtomSet iterates) {
        LinkedList list = generalIteratorMethodTests(treeIterator);
        Lister testLister = new Lister();
        testLister.addEachToList(iterates);
        assertEquals(list, testLister.list);
        return list;
    }
    
    private Phase phase;
    private SpeciesAgent speciesAgent0, speciesAgent1, speciesAgent2;
    private AtomIteratorTreePhase treeIterator;
    int n0a, nAtoms, n1a, n2a;
    int[] nTree;

}