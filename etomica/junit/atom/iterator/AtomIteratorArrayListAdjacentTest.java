package etomica.junit.atom.iterator;

import java.util.LinkedList;

import etomica.atom.AtomArrayList;
import etomica.atom.IAtom;
import etomica.atom.ISpeciesAgent;
import etomica.atom.iterator.AtomIteratorArrayListAdjacent;
import etomica.atom.iterator.IteratorDirective;
import etomica.junit.UnitTestUtil;
import etomica.simulation.ISimulation;


/**
 * Unit test for AtomIteratorSequenceAdjacent
 * 
 * @author David Kofke
 *
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
        ISimulation sim = UnitTestUtil.makeStandardSpeciesTree(
                new int[] {nAtoms},2,new int[] {nAtoms},null,null);
        ISpeciesAgent agent = sim.getBoxs()[0].getAgent(sim.getSpeciesManager().getSpecies()[0]);
        AtomArrayList atomList = new AtomArrayList();
        atomList.addAll(agent.getChildList());

        //atom in middle of list
        IAtom atom = atomList.getAtom(5);
        IAtom dnAtom = atomList.getAtom(4);
        IAtom upAtom = atomList.getAtom(6);
        LinkedList list = null;
        Lister test = new Lister();
        
        iteratorUp.setAtom(atom);
        iteratorDn.setAtom(atom);
        iteratorBoth.setAtom(atom);
        
        iterator = iteratorUp;
        list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(), 1);
        test.clear();
        test.addEachToList(upAtom);
        assertEquals(list, test.list);

        iterator = iteratorDn;
        list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(), 1);
        test.clear();
        test.addEachToList(dnAtom);
        assertEquals(list, test.list);

        iterator = iteratorBoth;
        list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(), 2);
        test.clear();
        test.addEachToList(upAtom);
        test.addEachToList(dnAtom);
        assertEquals(list, test.list);

        //atom at end of list        
        atom = atomList.getAtom(nAtoms-1);
        dnAtom = atomList.getAtom(nAtoms-2);
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
        test.addEachToList(dnAtom);
        assertEquals(list, test.list);

        iterator = iteratorBoth;
        list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(), 1);
        test.clear();
        test.addEachToList(dnAtom);
        assertEquals(list, test.list);

        //atom at beginning of list        
        atom = atomList.getAtom(0);
        upAtom = atomList.getAtom(1);
        iteratorUp.setAtom(atom);
        iteratorDn.setAtom(atom);
        iteratorBoth.setAtom(atom);
        
        iterator = iteratorUp;
        list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(), 1);
        test.clear();
        test.addEachToList(upAtom);
        assertEquals(list, test.list);

        iterator = iteratorDn;
        list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(), 0);

        iterator = iteratorBoth;
        list = generalIteratorMethodTests(iterator);
        assertEquals(list.size(), 1);
        test.clear();
        test.addEachToList(upAtom);
        assertEquals(list, test.list);

        //short list
        agent.setNMolecules(1);
        atomList = new AtomArrayList();
        atomList.addAll(agent.getChildList());
        assertEquals(atomList.getAtomCount(), 1);
        atom = atomList.getAtom(0);
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
