package etomica.junit.atom.iterator;

import java.util.LinkedList;

import etomica.atom.Atom;
import etomica.atom.AtomArrayList;
import etomica.atom.SpeciesAgent;
import etomica.atom.SpeciesRoot;
import etomica.atom.iterator.AtomIteratorArrayListAdjacent;
import etomica.atom.iterator.IteratorDirective;
import etomica.junit.UnitTestUtil;


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
public class AtomIteratorArrayListAdjacentTest extends IteratorTestAbstract {

    public AtomIteratorArrayListAdjacentTest() {
        super();
        UnitTestUtil.VERBOSE = false;
    }

    public void setUp() {
        iteratorUp = new AtomIteratorArrayListAdjacent(IteratorDirective.Direction.UP);
        iteratorDn = new AtomIteratorArrayListAdjacent(IteratorDirective.Direction.DOWN);
        iteratorBoth = new AtomIteratorArrayListAdjacent(null);
    }
    
    public void testIteration() {
        int nAtoms = 11;
        SpeciesRoot root = UnitTestUtil.makeStandardSpeciesTree(
                new int[] {nAtoms},2,new int[] {nAtoms},null,null);
        SpeciesAgent agent = (SpeciesAgent)root.getDescendant(new int[] {0,0});
        AtomArrayList atomList = new AtomArrayList();
        atomList.addAll(agent.getChildList());

        //atom in middle of list
        Atom atom = atomList.get(5);
        Atom dnAtom = atomList.get(4);
        Atom upAtom = atomList.get(6);
        LinkedList list = null;
        Lister test = new Lister();
        
        iteratorUp.setAtom(atom);
        iteratorDn.setAtom(atom);
        iteratorBoth.setAtom(atom);
        
        iterator = iteratorUp;
        list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(), 1);
        test.clear();
        test.addEachToList(new Atom[] {upAtom});
        assertEquals(list, test.list);

        iterator = iteratorDn;
        list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(), 1);
        test.clear();
        test.addEachToList(new Atom[] {dnAtom});
        assertEquals(list, test.list);

        iterator = iteratorBoth;
        list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(), 2);
        test.clear();
        test.addEachToList(new Atom[] {upAtom,dnAtom});
        assertEquals(list, test.list);

        //atom at end of list        
        atom = atomList.get(nAtoms-1);
        dnAtom = atomList.get(nAtoms-2);
        iteratorUp.setAtom(atom);
        iteratorDn.setAtom(atom);
        iteratorBoth.setAtom(atom);
        
        iterator = iteratorUp;
        list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(), 0);

        iterator = iteratorDn;
        list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(), 1);
        test.clear();
        test.addEachToList(new Atom[] {dnAtom});
        assertEquals(list, test.list);

        iterator = iteratorBoth;
        list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(), 1);
        test.clear();
        test.addEachToList(new Atom[] {dnAtom});
        assertEquals(list, test.list);

        //atom at beginning of list        
        atom = atomList.get(0);
        upAtom = atomList.get(1);
        iteratorUp.setAtom(atom);
        iteratorDn.setAtom(atom);
        iteratorBoth.setAtom(atom);
        
        iterator = iteratorUp;
        list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(), 1);
        test.clear();
        test.addEachToList(new Atom[] {upAtom});
        assertEquals(list, test.list);

        iterator = iteratorDn;
        list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(), 0);

        iterator = iteratorBoth;
        list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(), 1);
        test.clear();
        test.addEachToList(new Atom[] {upAtom});
        assertEquals(list, test.list);

        //short list
        agent.setNMolecules(1);
        atomList = new AtomArrayList();
        atomList.addAll(agent.getChildList());
        assertEquals(atomList.size(), 1);
        atom = atomList.get(0);
        iteratorUp.setAtom(atom);
        iteratorDn.setAtom(atom);
        iteratorBoth.setAtom(atom);
        iterator = iteratorUp;
        list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(), 0);
        iterator = iteratorDn;
        list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(), 0);
        iterator = iteratorBoth;
        list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(), 0);
        
    }
    
    AtomIteratorArrayListAdjacent iteratorUp, iteratorDn, iteratorBoth, iterator;
}
