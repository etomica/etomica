package etomica.junit.atom.iterator;

import java.util.LinkedList;

import etomica.Atom;
import etomica.AtomSet;
import etomica.AtomTreeNodeGroup;
import etomica.SpeciesMaster;
import etomica.SpeciesRoot;
import etomica.atom.AtomList;
import etomica.atom.iterator.AtomIteratorTree;
import etomica.junit.UnitTest;

/**
 * Unit test for AtomIteratorTree
 */
public class AtomIteratorTreeTest extends IteratorTest {

    protected void setUp() {
        n0a = 1;
        nAtoms = 3;
        n1a = 0;
        n2a = 0;
        nTree = new int[] { 5, 4, 3 };
        root = UnitTest.makeStandardSpeciesTree(new int[] { n0a },
                nAtoms, new int[] { n1a }, new int[] { n2a }, nTree);
        rootNode = (AtomTreeNodeGroup) root.node;
        speciesMaster = (SpeciesMaster)rootNode.getDescendant(new int[] {0});

        iterator = new AtomIteratorTree();
    }
    
    //species 0: 5 molecules, each a group of three atoms
    //species 1: 10 no-group single-atom molecules
    //species 2: 3 molecules, each with 5 subgroups of 4 groups of 3 atoms
    public void testIterator() {
        Lister testLister = new Lister();
        LinkedList list = null;
        Atom iterate = null;
        AtomList iterates = null;
        Atom iterationRoot = null;
        int count = 0;

        //test initial iterator provides no iterates
        list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(), 0);
        
        //test iteration over all nodes
        count = 1 + 1 + 3 + n0a*(1 + nAtoms) + n1a*(1) 
                + n2a*(1 + nTree[0]*(1 + nTree[1]*(1 + nTree[2])));
        iterationRoot = root;
        iterator.setDoAllNodes(true);
        list = testIterateCount(iterationRoot, count);

        if (n0a>0) {
            iterationRoot = rootNode.getDescendant(new int[]{0,0,0,0});
            testOneIterate(iterationRoot,iterationRoot);
        }
        
        //test iteration of leaf atoms
        iterationRoot = root;
        iterator.setAsLeafIterator();
        count = n0a*nAtoms + n1a + n2a*nTree[0]*nTree[1]*nTree[2];
        list = testIterateCount(root, count);
        list = testListIterates(root, speciesMaster.atomList);
        
        testNoIterates(null);
    }
    
    private LinkedList testIterateCount(Atom root, int count) {
        iterator.setRoot(root);
        LinkedList list = generalIteratorMethodTests(iterator);
        System.out.println(list.size()+ " "+count);
        assertEquals(list.size(), count);
        return list;
    }
    
    private LinkedList testOneIterate(Atom root, Atom iterate) {
        iterator.setRoot(root);
        LinkedList list = generalIteratorMethodTests(iterator);
        Lister testLister = new Lister();
        testLister.actionPerformed(iterate);
        assertEquals(list, testLister.list);
        return list;
    }
    
    private LinkedList testListIterates(Atom root, AtomList iterates) {
        iterator.setRoot(root);
        LinkedList list = generalIteratorMethodTests(iterator);
        Lister testLister = new Lister();
        testLister.addEachToList(iterates);
        assertEquals(list, testLister.list);
        return list;
    }
    
    private void testNoIterates(Atom root) {
        iterator.setRoot(root);
        LinkedList list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(), 0);
    }


    private SpeciesRoot root;
    private SpeciesMaster speciesMaster;
    private AtomIteratorTree iterator;
    private AtomTreeNodeGroup rootNode;
    int n0a, nAtoms, n1a, n2a;
    int[] nTree;

}