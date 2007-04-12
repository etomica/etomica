package etomica.junit.atom.iterator;

import java.util.LinkedList;

import etomica.atom.AtomArrayList;
import etomica.atom.IAtom;
import etomica.atom.SpeciesAgent;
import etomica.atom.iterator.AtomIteratorTreeRoot;
import etomica.junit.UnitTestUtil;
import etomica.phase.Phase;
import etomica.simulation.Simulation;
import etomica.species.Species;

/**
 * Unit test for AtomIteratorTree
 */
public class AtomIteratorTreeTest extends IteratorTestAbstract {

    protected void setUp() {
        n0a = 3;
        nAtoms = 3;
        n1a = 10;
        n2a = 4;
        nTree = new int[] { 5, 4, 3 };
        Simulation sim = UnitTestUtil.makeStandardSpeciesTree(new int[] { n0a },
                nAtoms, new int[] { n1a }, new int[] { n2a }, nTree);
        Phase phase = sim.getPhases()[0];
        Species[] species = sim.getSpeciesManager().getSpecies();
        speciesAgent0 = phase.getAgent(species[0]);
        speciesAgent2 = phase.getAgent(species[2]);

        treeIterator = new AtomIteratorTreeRoot();
    }
    
    //species 0: 5 molecules, each a group of three atoms
    //species 1: 10 no-group single-atom molecules
    //species 2: 3 molecules, each with 5 subgroups of 4 groups of 3 atoms
    public void testIterator() {
        LinkedList list = null;
        IAtom iterationRoot = null;
        int count = 0;

        //test initial iterator provides no iterates
        list = generalIteratorMethodTests(treeIterator);
        assertEquals(list.size(), 0);
        
        //test iteration over all nodes from speciesAgent2
        count = 1 + n2a*(1 + nTree[0]*(1 + nTree[1]*(1 + nTree[2])));
        iterationRoot = speciesAgent2;
        treeIterator.setDoAllNodes(true);
        list = testIterateCount(iterationRoot, count);

        //test iteration over different depths, starting at root
        iterationRoot = speciesAgent2;
        treeIterator.setDoAllNodes(true);
        list = testOneIterate(0, iterationRoot, iterationRoot);
        testNoIterates((IAtom)null);
        list = testArrayIterates(0, iterationRoot, new IAtom[] {speciesAgent2});
        list = testIterateCount(1, iterationRoot, 1+n2a);
        list = testIterateCount(2, iterationRoot, 1+n2a*(1+nTree[0]));
        list = testIterateCount(3, iterationRoot, 1+n2a*(1+nTree[0]*(1+nTree[1])));
        list = testIterateCount(4, iterationRoot, 1+n2a*(1+nTree[0]*(1+nTree[1]*(1+nTree[2]))));
        list = testIterateCount(100, iterationRoot, 1+n2a*(1+nTree[0]*(1+nTree[1]*(1+nTree[2]))));
        testNoIterates((IAtom)null);

        treeIterator.setDoAllNodes(false);
        list = testOneIterate(0, iterationRoot, iterationRoot);
        testNoIterates((IAtom)null);
        list = testArrayIterates(0, iterationRoot, new IAtom[] {speciesAgent2});
        list = testIterateCount(1, iterationRoot, n2a);
        list = testIterateCount(2, iterationRoot, n2a*nTree[0]);
        list = testIterateCount(3, iterationRoot, n2a*nTree[0]*(nTree[1]));
        list = testIterateCount(4, iterationRoot, n2a*nTree[0]*nTree[1]*nTree[2]);
        list = testIterateCount(100, iterationRoot, n2a*nTree[0]*nTree[1]*nTree[2]);
        testNoIterates((IAtom)null);

        //test iteration over different depths, starting at a species agent
        iterationRoot = speciesAgent2;
        treeIterator.setDoAllNodes(true);
        list = testOneIterate(0, iterationRoot, iterationRoot);
        testNoIterates((IAtom)null);
        AtomArrayList testList = (AtomArrayList)speciesAgent2.getChildList().clone();
        testList.add(0,iterationRoot);
        list = testListIterates(1, iterationRoot, testList);
        list = testIterateCount(2, iterationRoot, 1+n2a*(1+nTree[0]));
        list = testIterateCount(3, iterationRoot, 1+n2a*(1+nTree[0]*(1+nTree[1])));
        list = testIterateCount(4, iterationRoot, 1+n2a*(1+nTree[0]*(1+nTree[1]*(1+nTree[2]))));
        list = testIterateCount(100, iterationRoot, 1+n2a*(1+nTree[0]*(1+nTree[1]*(1+nTree[2]))));
        testNoIterates((IAtom)null);

        treeIterator.setDoAllNodes(false);
        list = testOneIterate(0, iterationRoot, iterationRoot);
        testNoIterates((IAtom)null);
        testList = (AtomArrayList)speciesAgent2.getChildList().clone();
        list = testListIterates(1, iterationRoot, testList);
        list = testIterateCount(2, iterationRoot, n2a*nTree[0]);
        list = testIterateCount(3, iterationRoot, n2a*nTree[0]*(nTree[1]));
        list = testIterateCount(4, iterationRoot, n2a*nTree[0]*nTree[1]*nTree[2]);
        list = testIterateCount(100, iterationRoot, n2a*nTree[0]*nTree[1]*nTree[2]);
        testNoIterates((IAtom)null);
        
        //test leaf atom selected as root
        if (n0a>0) {
            iterationRoot = speciesAgent0.getDescendant(new int[]{0,0});
            treeIterator.setDoAllNodes(true);
            testOneIterate(iterationRoot,iterationRoot);
            treeIterator.setDoAllNodes(false);
            testOneIterate(iterationRoot,iterationRoot);
        }
        
        //test iteration of leaf atoms
        iterationRoot = speciesAgent2;
        treeIterator.setAsLeafIterator();
        count = n2a*nTree[0]*nTree[1]*nTree[2];
        list = testIterateCount(iterationRoot, count);
        
        //test re-specifying iteration in different orders
        iterationRoot = speciesAgent2;
        int depth = 2;
        boolean doAllNodes = true;
        count = 1+n2a*(1+nTree[0]);
        treeIterator = new AtomIteratorTreeRoot(speciesAgent2, depth, doAllNodes);
        list = generalIteratorMethodTests(treeIterator);
        assertEquals(list.size(), count);
        treeIterator = new AtomIteratorTreeRoot();
        treeIterator.setDoAllNodes(doAllNodes);
        treeIterator.setRootAtom(speciesAgent2);
        treeIterator.setIterationDepth(depth);
        LinkedList list1 = generalIteratorMethodTests(treeIterator);
        assertEquals(list, list1);

        treeIterator.setIterationDepth(0);
        treeIterator.setRootAtom(null);
        treeIterator.setDoAllNodes(!doAllNodes);
        treeIterator.setIterationDepth(depth);
        treeIterator.setDoAllNodes(doAllNodes);
        treeIterator.setRootAtom(speciesAgent2);
        list1 = generalIteratorMethodTests(treeIterator);
        assertEquals(list, list1);
        
        treeIterator.setIterationDepth(0);
        treeIterator.setRootAtom(null);
        treeIterator.setDoAllNodes(!doAllNodes);
        treeIterator.setRootAtom(speciesAgent2);
        treeIterator.setIterationDepth(depth);
        treeIterator.setDoAllNodes(doAllNodes);
        list1 = generalIteratorMethodTests(treeIterator);
        assertEquals(list, list1);

        
        testNoIterates((IAtom)null);
    }
    
    private LinkedList testIterateCount(int depth, IAtom iterationRoot, int count) {
        treeIterator.setIterationDepth(depth);
        LinkedList list0 = testIterateCount(iterationRoot, count);
        treeIterator.setIterationDepth(0);
        treeIterator.setRootAtom(null);
        treeIterator.setRootAtom(iterationRoot);//change order of setting
        treeIterator.setIterationDepth(depth);
        LinkedList list1 = generalIteratorMethodTests(treeIterator);
        assertEquals(list0, list1);
        return list0;
    }
    private LinkedList testIterateCount(IAtom iterationRoot, int count) {
        treeIterator.setRootAtom(iterationRoot);
        LinkedList list = generalIteratorMethodTests(treeIterator);
//        System.out.println(list.size()+ " "+count);
        assertEquals(list.size(), count);
        return list;
    }
    
    private LinkedList testOneIterate(int depth, IAtom iterationRoot, IAtom iterate) {
        treeIterator.setIterationDepth(depth);
        return testOneIterate(iterationRoot, iterate);
    }
    private LinkedList testOneIterate(IAtom iterationRoot, IAtom iterate) {
        treeIterator.setRootAtom(iterationRoot);
        LinkedList list = generalIteratorMethodTests(treeIterator);
        Lister testLister = new Lister();
        testLister.actionPerformed(iterate);
        assertEquals(list, testLister.list);
        return list;
    }
    
    private LinkedList testArrayIterates(int depth, IAtom iterationRoot, IAtom[] iterates) {
        treeIterator.setIterationDepth(depth);
        AtomArrayList arrayList = new AtomArrayList();
        for (int i=0; i<iterates.length; i++) {
            arrayList.add(iterates[i]);
        }
        return testListIterates(iterationRoot, arrayList);
    }
    
    private LinkedList testListIterates(int depth, IAtom iterationRoot, AtomArrayList iterates) {
        treeIterator.setIterationDepth(depth);
        return testListIterates(iterationRoot, iterates);
    }
    private LinkedList testListIterates(IAtom iterationRoot, AtomArrayList iterates) {
        treeIterator.setRootAtom(iterationRoot);
        LinkedList list = generalIteratorMethodTests(treeIterator);
        Lister testLister = new Lister();
        testLister.addEachToList(iterates);
        assertEquals(list, testLister.list);
        return list;
    }
    
    private void testNoIterates(IAtom iterationRoot) {
        treeIterator.setRootAtom(iterationRoot);
        LinkedList list = generalIteratorMethodTests(treeIterator);
        assertEquals(list.size(), 0);
    }

    private SpeciesAgent speciesAgent0, speciesAgent2;
    private AtomIteratorTreeRoot treeIterator;
    int n0a, nAtoms, n1a, n2a;
    int[] nTree;

}