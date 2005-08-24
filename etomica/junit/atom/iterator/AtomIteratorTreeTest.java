package etomica.junit.atom.iterator;

import java.util.LinkedList;

import etomica.atom.Atom;
import etomica.atom.AtomList;
import etomica.atom.AtomTreeNodeGroup;
import etomica.atom.SpeciesAgent;
import etomica.atom.SpeciesMaster;
import etomica.atom.SpeciesRoot;
import etomica.atom.iterator.AtomIteratorTree;
import etomica.junit.UnitTest;

/**
 * Unit test for AtomIteratorTree
 */
public class AtomIteratorTreeTest extends IteratorTest {

    protected void setUp() {
        n0a = 3;
        nAtoms = 3;
        n1a = 10;
        n2a = 4;
        nTree = new int[] { 5, 4, 3 };
        root = UnitTest.makeStandardSpeciesTree(new int[] { n0a },
                nAtoms, new int[] { n1a }, new int[] { n2a }, nTree);
        rootNode = (AtomTreeNodeGroup) root.node;
        speciesMaster = (SpeciesMaster)rootNode.getDescendant(new int[] {0});
        speciesAgent0 = (SpeciesAgent)rootNode.getDescendant(new int[] {0,0});
        speciesAgent1 = (SpeciesAgent)rootNode.getDescendant(new int[] {0,1});
        speciesAgent2 = (SpeciesAgent)rootNode.getDescendant(new int[] {0,2});

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
        
        //test iteration over all nodes from root
        count = 1 + 1 + 3 + n0a*(1 + nAtoms) + n1a*(1) 
                + n2a*(1 + nTree[0]*(1 + nTree[1]*(1 + nTree[2])));
        iterationRoot = root;
        treeIterator.setDoAllNodes(true);
        list = testIterateCount(iterationRoot, count);

        //test iteration over different depths, starting at root
        iterationRoot = root;
        treeIterator.setDoAllNodes(true);
        list = testOneIterate(0, iterationRoot, iterationRoot);
        testNoIterates((Atom)null);
        list = testArrayIterates(1, iterationRoot, new Atom[] {iterationRoot, speciesMaster});
        list = testArrayIterates(2, iterationRoot, new Atom[] {iterationRoot, speciesMaster, speciesAgent0, speciesAgent1, speciesAgent2});
        list = testIterateCount(3, iterationRoot, 1+1+3+n0a+n1a+n2a);
        list = testIterateCount(4, iterationRoot, 1+1+3+n0a*(1+nAtoms)+n1a+n2a*(1+nTree[0]));
        list = testIterateCount(5, iterationRoot, 1+1+3+n0a*(1+nAtoms)+n1a+n2a*(1+nTree[0]*(1+nTree[1])));
        list = testIterateCount(6, iterationRoot, 1+1+3+n0a*(1+nAtoms)+n1a+n2a*(1+nTree[0]*(1+nTree[1]*(1+nTree[2]))));
        list = testIterateCount(100, iterationRoot, 1+1+3+n0a*(1+nAtoms)+n1a+n2a*(1+nTree[0]*(1+nTree[1]*(1+nTree[2]))));
        testNoIterates((Atom)null);

        treeIterator.setDoAllNodes(false);
        list = testOneIterate(0, iterationRoot, iterationRoot);
        testNoIterates((Atom)null);
        list = testOneIterate(1, iterationRoot, speciesMaster);
        list = testArrayIterates(2, iterationRoot, new Atom[] {speciesAgent0, speciesAgent1, speciesAgent2});
        list = testIterateCount(3, iterationRoot, n0a+n1a+n2a);
        list = testIterateCount(4, iterationRoot, n0a*nAtoms+n1a+n2a*nTree[0]);
        list = testIterateCount(5, iterationRoot, n0a*nAtoms+n1a+n2a*nTree[0]*(nTree[1]));
        list = testIterateCount(6, iterationRoot, n0a*nAtoms+n1a+n2a*nTree[0]*nTree[1]*nTree[2]);
        list = testIterateCount(100, iterationRoot, n0a*nAtoms+n1a+n2a*nTree[0]*nTree[1]*nTree[2]);
        testNoIterates((Atom)null);

        //test iteration over different depths, starting at a species agent
        iterationRoot = speciesAgent2;
        treeIterator.setDoAllNodes(true);
        list = testOneIterate(0, iterationRoot, iterationRoot);
        testNoIterates((Atom)null);
        AtomList testList = new AtomList(((AtomTreeNodeGroup)speciesAgent2.node).childList);
        testList.addFirst(iterationRoot);
        list = testListIterates(1, iterationRoot, testList);
        list = testIterateCount(2, iterationRoot, 1+n2a*(1+nTree[0]));
        list = testIterateCount(3, iterationRoot, 1+n2a*(1+nTree[0]*(1+nTree[1])));
        list = testIterateCount(4, iterationRoot, 1+n2a*(1+nTree[0]*(1+nTree[1]*(1+nTree[2]))));
        list = testIterateCount(100, iterationRoot, 1+n2a*(1+nTree[0]*(1+nTree[1]*(1+nTree[2]))));
        testNoIterates((Atom)null);

        treeIterator.setDoAllNodes(false);
        list = testOneIterate(0, iterationRoot, iterationRoot);
        testNoIterates((Atom)null);
        testList = ((AtomTreeNodeGroup)speciesAgent2.node).childList;
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
            iterationRoot = rootNode.getDescendant(new int[]{0,0,0,0});
            treeIterator.setDoAllNodes(true);
            testOneIterate(iterationRoot,iterationRoot);
            treeIterator.setDoAllNodes(false);
            testOneIterate(iterationRoot,iterationRoot);
        }
        
        //test iteration of leaf atoms
        iterationRoot = root;
        treeIterator.setAsLeafIterator();
        count = n0a*nAtoms + n1a + n2a*nTree[0]*nTree[1]*nTree[2];
        list = testIterateCount(root, count);
        list = testListIterates(root, speciesMaster.atomList);
        
        //test re-specifying iteration in different orders
        iterationRoot = root;
        int depth = 4;
        boolean doAllNodes = true;
        count = 1+1+3+n0a*(1+nAtoms)+n1a+n2a*(1+nTree[0]);
        treeIterator = new AtomIteratorTree(root, depth, doAllNodes);
        list = generalIteratorMethodTests(treeIterator);
        assertEquals(list.size(), count);
        treeIterator = new AtomIteratorTree();
        treeIterator.setDoAllNodes(doAllNodes);
        treeIterator.setRoot(root);
        treeIterator.setIterationDepth(depth);
        LinkedList list1 = generalIteratorMethodTests(treeIterator);
        assertEquals(list, list1);

        treeIterator.setIterationDepth(0);
        treeIterator.setRoot(null);
        treeIterator.setDoAllNodes(!doAllNodes);
        treeIterator.setIterationDepth(depth);
        treeIterator.setDoAllNodes(doAllNodes);
        treeIterator.setRoot(root);
        list1 = generalIteratorMethodTests(treeIterator);
        assertEquals(list, list1);
        
        treeIterator.setIterationDepth(0);
        treeIterator.setRoot(null);
        treeIterator.setDoAllNodes(!doAllNodes);
        treeIterator.setRoot(root);
        treeIterator.setIterationDepth(depth);
        treeIterator.setDoAllNodes(doAllNodes);
        list1 = generalIteratorMethodTests(treeIterator);
        assertEquals(list, list1);

        
        testNoIterates((Atom)null);
    }
    
    private LinkedList testIterateCount(int depth, Atom root, int count) {
        treeIterator.setIterationDepth(depth);
        LinkedList list0 = testIterateCount(root, count);
        treeIterator.setIterationDepth(0);
        treeIterator.setRoot(null);
        treeIterator.setRoot(root);//change order of setting
        treeIterator.setIterationDepth(depth);
        LinkedList list1 = generalIteratorMethodTests(treeIterator);
        assertEquals(list0, list1);
        return list0;
    }
    private LinkedList testIterateCount(Atom root, int count) {
        treeIterator.setRoot(root);
        LinkedList list = generalIteratorMethodTests(treeIterator);
//        System.out.println(list.size()+ " "+count);
        assertEquals(list.size(), count);
        return list;
    }
    
    private LinkedList testOneIterate(int depth, Atom root, Atom iterate) {
        treeIterator.setIterationDepth(depth);
        return testOneIterate(root, iterate);
    }
    private LinkedList testOneIterate(Atom root, Atom iterate) {
        treeIterator.setRoot(root);
        LinkedList list = generalIteratorMethodTests(treeIterator);
        Lister testLister = new Lister();
        testLister.actionPerformed(iterate);
        assertEquals(list, testLister.list);
        return list;
    }
    
    private LinkedList testArrayIterates(int depth, Atom root, Atom[] iterates) {
        treeIterator.setIterationDepth(depth);
        return testListIterates(root, new AtomList(iterates));
    }
    
    private LinkedList testListIterates(int depth, Atom root, AtomList iterates) {
        treeIterator.setIterationDepth(depth);
        return testListIterates(root, iterates);
    }
    private LinkedList testListIterates(Atom root, AtomList iterates) {
        treeIterator.setRoot(root);
        LinkedList list = generalIteratorMethodTests(treeIterator);
        Lister testLister = new Lister();
        testLister.addEachToList(iterates);
        assertEquals(list, testLister.list);
        return list;
    }
    
    private void testNoIterates(Atom root) {
        treeIterator.setRoot(root);
        LinkedList list = generalIteratorMethodTests(treeIterator);
        assertEquals(list.size(), 0);
    }


    private SpeciesRoot root;
    private SpeciesMaster speciesMaster;
    private SpeciesAgent speciesAgent0, speciesAgent1, speciesAgent2;
    private AtomIteratorTree treeIterator;
    private AtomTreeNodeGroup rootNode;
    int n0a, nAtoms, n1a, n2a;
    int[] nTree;

}