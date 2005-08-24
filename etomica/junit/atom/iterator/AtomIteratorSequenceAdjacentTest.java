package etomica.junit.atom.iterator;

import java.util.LinkedList;

import etomica.IteratorDirective;
import etomica.atom.Atom;
import etomica.atom.AtomList;
import etomica.atom.AtomTreeNodeGroup;
import etomica.atom.SpeciesAgent;
import etomica.atom.SpeciesRoot;
import etomica.atom.iterator.AtomIteratorSequenceAdjacent;
import etomica.junit.UnitTest;


/**
 * Unit test for AtomIteratorSequenceAdjacent
 * 
 * @author David Kofke
 *
 */

/*
 * History
 * Created on Jun 14, 2005 by kofke
 */
public class AtomIteratorSequenceAdjacentTest extends IteratorTest {

    public AtomIteratorSequenceAdjacentTest() {
        super();
        UnitTest.VERBOSE = false;
    }

    public void setUp() {
        iterator = new AtomIteratorSequenceAdjacent();
    }
    
    public void testIteration() {
        int nAtoms = 11;
        SpeciesRoot root = UnitTest.makeStandardSpeciesTree(
                new int[] {nAtoms},2,new int[] {nAtoms},null,null);
        AtomTreeNodeGroup rootNode = (AtomTreeNodeGroup)root.node;
        SpeciesAgent agent = (SpeciesAgent)rootNode.getDescendant(new int[] {0,0});
        AtomList atomList = ((AtomTreeNodeGroup)agent.node).childList;

        //atom in middle of list
        Atom atom = atomList.get(5);
        Atom dnAtom = atomList.get(4);
        Atom upAtom = atomList.get(6);
        iterator.setAtom(atom);
        iterator.setDirection(null);
        LinkedList list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(), 2);
        Lister test = new Lister();
        test.addEachToList(new Atom[] {upAtom, dnAtom});
        assertEquals(list, test.list);
        
        iterator.setDirection(IteratorDirective.UP);
        list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(), 1);
        test.clear();
        test.addEachToList(new Atom[] {upAtom});
        assertEquals(list, test.list);

        iterator.setDirection(IteratorDirective.DOWN);
        list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(), 1);
        test.clear();
        test.addEachToList(new Atom[] {dnAtom});
        assertEquals(list, test.list);

        //atom at end of list        
        atom = atomList.get(nAtoms-1);
        dnAtom = atomList.get(nAtoms-2);
        iterator.setAtom(atom);
        iterator.setDirection(null);
        list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(), 1);
        test.clear();
        test.addEachToList(new Atom[] {dnAtom});
        assertEquals(list, test.list);
        
        iterator.setDirection(IteratorDirective.UP);
        list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(), 0);

        iterator.setDirection(IteratorDirective.DOWN);
        list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(), 1);
        test.clear();
        test.addEachToList(new Atom[] {dnAtom});
        assertEquals(list, test.list);

        //atom at beginning of list        
        atom = atomList.get(0);
        upAtom = atomList.get(1);
        iterator.setAtom(atom);
        iterator.setDirection(null);
        list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(), 1);
        test.clear();
        test.addEachToList(new Atom[] {upAtom});
        assertEquals(list, test.list);
        
        iterator.setDirection(IteratorDirective.UP);
        list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(), 1);
        test.clear();
        test.addEachToList(new Atom[] {upAtom});
        assertEquals(list, test.list);

        iterator.setDirection(IteratorDirective.DOWN);
        list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(), 0);

        //short list
        agent.setNMolecules(1);
        assertEquals(atomList.size(), 1);
        atom = atomList.get(0);
        iterator.setAtom(atom);
        iterator.setDirection(null);
        list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(), 0);
        iterator.setDirection(IteratorDirective.UP);
        list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(), 0);
        iterator.setDirection(IteratorDirective.DOWN);
        list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(), 0);
        
    }
    
    AtomIteratorSequenceAdjacent iterator;
}
