package etomica.junit.atom.iterator;

import java.util.LinkedList;

import etomica.atom.Atom;
import etomica.atom.AtomArrayList;
import etomica.atom.SpeciesAgent;
import etomica.atom.SpeciesMaster;
import etomica.atom.iterator.AtomIteratorTree;
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
        speciesMaster = phase.getSpeciesMaster();
        Species[] species = sim.getSpeciesManager().getSpecies();
        speciesAgent0 = phase.getAgent(species[0]);
        speciesAgent1 = phase.getAgent(species[1]);
        speciesAgent2 = phase.getAgent(species[2]);

        treeIterator = new AtomIteratorTree();
    }
    
    //species 0: 5 molecules, each a group of three atoms
    //species 1: 10 no-group single-atom molecules
    //species 2: 3 molecules, each with 5 subgroups of 4 groups of 3 atoms
    public void testIterator() {
        LinkedList list = null;
        Atom iterationRoot = null;
        int count = 0;

        //test initial iterator provides no iterates
        list = generalIteratorMethodTests(treeIterator);
        assertEquals(list.size(), 0);
        
        //test iteration over all nodes from speciesMaster
        count = 1 + 3 + n0a*(1 + nAtoms) + n1a*(1) 
                + n2a*(1 + nTree[0]*(1 + nTree[1]*(1 + nTree[2])));
        iterationRoot = speciesMaster;
        treeIterator.setDoAllNodes(true);
        list = testIterateCount(iterationRoot, count);

        //test iteration over different depths, starting at root
        iterationRoot = speciesMaster;
        treeIterator.setDoAllNodes(true);
        list = testOneIterate(0, iterationRoot, iterationRoot);
        testNoIterates((Atom)null);
        list = testArrayIterates(1, iterationRoot, new Atom[] {speciesMaster, speciesAgent0, speciesAgent1, speciesAgent2});
        list = testIterateCount(2, iterationRoot, 1+3+n0a+n1a+n2a);
        list = testIterateCount(3, iterationRoot, 1+3+n0a*(1+nAtoms)+n1a+n2a*(1+nTree[0]));
        list = testIterateCount(4, iterationRoot, 1+3+n0a*(1+nAtoms)+n1a+n2a*(1+nTree[0]*(1+nTree[1])));
        list = testIterateCount(5, iterationRoot, 1+3+n0a*(1+nAtoms)+n1a+n2a*(1+nTree[0]*(1+nTree[1]*(1+nTree[2]))));
        list = testIterateCount(100, iterationRoot, 1+3+n0a*(1+nAtoms)+n1a+n2a*(1+nTree[0]*(1+nTree[1]*(1+nTree[2]))));
        testNoIterates((Atom)null);

        treeIterator.setDoAllNodes(false);
        list = testOneIterate(0, iterationRoot, iterationRoot);
        testNoIterates((Atom)null);
        list = testArrayIterates(1, iterationRoot, new Atom[] {speciesAgent0, speciesAgent1, speciesAgent2});
        list = testIterateCount(2, iterationRoot, n0a+n1a+n2a);
        list = testIterateCount(3, iterationRoot, n0a*nAtoms+n1a+n2a*nTree[0]);
        list = testIterateCount(4, iterationRoot, n0a*nAtoms+n1a+n2a*nTree[0]*(nTree[1]));
        list = testIterateCount(5, iterationRoot, n0a*nAtoms+n1a+n2a*nTree[0]*nTree[1]*nTree[2]);
        list = testIterateCount(100, iterationRoot, n0a*nAtoms+n1a+n2a*nTree[0]*nTree[1]*nTree[2]);
        testNoIterates((Atom)null);

        //test iteration over different depths, starting at a species agent
        iterationRoot = speciesAgent2;
        treeIterator.setDoAllNodes(true);
        list = testOneIterate(0, iterationRoot, iterationRoot);
        testNoIterates((Atom)null);
        AtomArrayList testList = (AtomArrayList)speciesAgent2.getChildList().clone();
        testList.add(0,iterationRoot);
        list = testListIterates(1, iterationRoot, testList);
        list = testIterateCount(2, iterationRoot, 1+n2a*(1+nTree[0]));
        list = testIterateCount(3, iterationRoot, 1+n2a*(1+nTree[0]*(1+nTree[1])));
        list = testIterateCount(4, iterationRoot, 1+n2a*(1+nTree[0]*(1+nTree[1]*(1+nTree[2]))));
        list = testIterateCount(100, iterationRoot, 1+n2a*(1+nTree[0]*(1+nTree[1]*(1+nTree[2]))));
        testNoIterates((Atom)null);

        treeIterator.setDoAllNodes(false);
        list = testOneIterate(0, iterationRoot, iterationRoot);
        testNoIterates((Atom)null);
        testList = (AtomArrayList)speciesAgent2.getChildList().clone();
        list = testListIterates(1, iterationRoot, testList);
        list = testIterateCount(2, iterationRoot, n2a*nTree[0]);
        list = testIterateCount(3, iterationRoot, n2a*nTree[0]*(nTree[1]));
        list = testIterateCount(4, iterationRoot, n2a*nTree[0]*nTree[1]*nTree[2]);
        list = testIterateCount(100, iterationRoot, n2a*nTree[0]*nTree[1]*nTree[2]);
        testNoIterates((Atom)null);
        
        count = 1 + 1 + 3 + n0a*(1 + nAtoms) + n1a*(1) 
        + n2a*(1 + nTree[0]*(1 + nTree[1]*(1 + nTree[2])));
        
        //test leaf atom selected as root
        if (n0a>0) {
            iterationRoot = speciesMaster.getDescendant(new int[]{0,0,0});
            treeIterator.setDoAllNodes(true);
            testOneIterate(iterationRoot,iterationRoot);
            treeIterator.setDoAllNodes(false);
            testOneIterate(iterationRoot,iterationRoot);
        }
        
        //test iteration of leaf atoms
        iterationRoot = speciesMaster;
        treeIterator.setAsLeafIterator();
        count = n0a*nAtoms + n1a + n2a*nTree[0]*nTree[1]*nTree[2];
        list = testIterateCount(speciesMaster, count);
        list = testListIterates(speciesMaster, speciesMaster.leafList);
        
        //test re-specifying iteration in different orders
        iterationRoot = speciesMaster;
        int depth = 3;
        boolean doAllNodes = true;
        count = 1+3+n0a*(1+nAtoms)+n1a+n2a*(1+nTree[0]);
        treeIterator = new AtomIteratorTree(speciesMaster, depth, doAllNodes);
        list = generalIteratorMethodTests(treeIterator);
        assertEquals(list.size(), count);
        treeIterator = new AtomIteratorTree();
        treeIterator.setDoAllNodes(doAllNodes);
        treeIterator.setRootAtom(speciesMaster);
        treeIterator.setIterationDepth(depth);
        LinkedList list1 = generalIteratorMethodTests(treeIterator);
        assertEquals(list, list1);

        treeIterator.setIterationDepth(0);
        treeIterator.setRootAtom(null);
        treeIterator.setDoAllNodes(!doAllNodes);
        treeIterator.setIterationDepth(depth);
        treeIterator.setDoAllNodes(doAllNodes);
        treeIterator.setRootAtom(speciesMaster);
        list1 = generalIteratorMethodTests(treeIterator);
        assertEquals(list, list1);
        
        treeIterator.setIterationDepth(0);
        treeIterator.setRootAtom(null);
        treeIterator.setDoAllNodes(!doAllNodes);
        treeIterator.setRootAtom(speciesMaster);
        treeIterator.setIterationDepth(depth);
        treeIterator.setDoAllNodes(doAllNodes);
        list1 = generalIteratorMethodTests(treeIterator);
        assertEquals(list, list1);

        
        testNoIterates((Atom)null);
    }
    
    private LinkedList testIterateCount(int depth, Atom iterationRoot, int count) {
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
    private LinkedList testIterateCount(Atom iterationRoot, int count) {
        treeIterator.setRootAtom(iterationRoot);
        LinkedList list = generalIteratorMethodTests(treeIterator);
//        System.out.println(list.size()+ " "+count);
        assertEquals(list.size(), count);
        return list;
    }
    
    private LinkedList testOneIterate(int depth, Atom iterationRoot, Atom iterate) {
        treeIterator.setIterationDepth(depth);
        return testOneIterate(iterationRoot, iterate);
    }
    private LinkedList testOneIterate(Atom iterationRoot, Atom iterate) {
        treeIterator.setRootAtom(iterationRoot);
        LinkedList list = generalIteratorMethodTests(treeIterator);
        Lister testLister = new Lister();
        testLister.actionPerformed(iterate);
        assertEquals(list, testLister.list);
        return list;
    }
    
    private LinkedList testArrayIterates(int depth, Atom iterationRoot, Atom[] iterates) {
        treeIterator.setIterationDepth(depth);
        AtomArrayList arrayList = new AtomArrayList();
        for (int i=0; i<iterates.length; i++) {
            arrayList.add(iterates[i]);
        }
        return testListIterates(iterationRoot, arrayList);
    }
    
    private LinkedList testListIterates(int depth, Atom iterationRoot, AtomArrayList iterates) {
        treeIterator.setIterationDepth(depth);
        return testListIterates(iterationRoot, iterates);
    }
    private LinkedList testListIterates(Atom iterationRoot, AtomArrayList iterates) {
        treeIterator.setRootAtom(iterationRoot);
        LinkedList list = generalIteratorMethodTests(treeIterator);
        Lister testLister = new Lister();
        testLister.addEachToList(iterates);
        assertEquals(list, testLister.list);
        return list;
    }
    
    private void testNoIterates(Atom iterationRoot) {
        treeIterator.setRootAtom(iterationRoot);
        LinkedList list = generalIteratorMethodTests(treeIterator);
        assertEquals(list.size(), 0);
    }

    private SpeciesMaster speciesMaster;
    private SpeciesAgent speciesAgent0, speciesAgent1, speciesAgent2;
    private AtomIteratorTree treeIterator;
    int n0a, nAtoms, n1a, n2a;
    int[] nTree;

}