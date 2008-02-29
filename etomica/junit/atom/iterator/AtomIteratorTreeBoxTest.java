package etomica.junit.atom.iterator;

import java.util.LinkedList;

import etomica.api.IBox;
import etomica.atom.AtomSet;
import etomica.atom.iterator.AtomIteratorTreeBox;
import etomica.junit.UnitTestUtil;
import etomica.simulation.ISimulation;

/**
 * Unit test for AtomIteratorTree
 */
public class AtomIteratorTreeBoxTest extends IteratorTestAbstract {

    protected void setUp() {
        n0a = 3;
        nAtoms = 3;
        n1a = 10;
        ISimulation sim = UnitTestUtil.makeStandardSpeciesTree(new int[] { n0a },
                nAtoms, new int[] { n1a });
        box = sim.getBoxs()[0];

        treeIterator = new AtomIteratorTreeBox();
        treeIterator.setBox(box);
    }
    
    //species 0: 5 molecules, each a group of three atoms
    //species 1: 10 no-group single-atom molecules
    //species 2: 3 molecules, each with 5 subgroups of 4 groups of 3 atoms
    public void testIterator() {
        LinkedList list = null;
        int count = 0;

        //test iteration over all nodes from speciesMaster
        count = n0a*(1 + nAtoms) + n1a*(1+1);
        treeIterator.setDoAllNodes(true);
        list = testIterateCount(count);

        //test iteration over different depths, starting at root
        treeIterator.setDoAllNodes(true);

        list = testIterateCount(1, n0a+n1a);
        list = testIterateCount(2, n0a*(1+nAtoms)+n1a*(1+1));

        treeIterator.setDoAllNodes(false);
        list = testIterateCount(1, n0a+n1a);
        list = testIterateCount(2, n0a*nAtoms+n1a);

        //test iteration of leaf atoms
        treeIterator.setAsLeafIterator();
        count = n0a*nAtoms + n1a;
        list = testIterateCount(count);
        list = testListIterates(box.getLeafList());
        
        //test re-specifying iteration in different orders
        int depth = 3;
        boolean doAllNodes = true;
        count = n0a*(1+nAtoms)+n1a*(1+1);
        treeIterator = new AtomIteratorTreeBox(box, depth, doAllNodes);
        list = generalIteratorMethodTests(treeIterator);
        assertEquals(list.size(), count);
        treeIterator = new AtomIteratorTreeBox();
        treeIterator.setDoAllNodes(doAllNodes);
        treeIterator.setBox(box);
        treeIterator.setIterationDepth(depth);
        LinkedList list1 = generalIteratorMethodTests(treeIterator);
        assertEquals(list, list1);

        treeIterator.setIterationDepth(depth);
        treeIterator.setDoAllNodes(doAllNodes);
        treeIterator.setBox(box);
        list1 = generalIteratorMethodTests(treeIterator);
        assertEquals(list, list1);
        
        treeIterator.setBox(box);
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
    
    private LinkedList testListIterates(AtomSet iterates) {
        LinkedList list = generalIteratorMethodTests(treeIterator);
        Lister testLister = new Lister();
        testLister.addEachToList(iterates);
        assertEquals(list, testLister.list);
        return list;
    }
    
    private IBox box;
    private AtomIteratorTreeBox treeIterator;
    int n0a, nAtoms, n1a;
}