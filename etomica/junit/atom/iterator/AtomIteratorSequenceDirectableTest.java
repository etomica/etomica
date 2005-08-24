package etomica.junit.atom.iterator;

import java.util.LinkedList;

import etomica.atom.AtomList;
import etomica.atom.AtomTreeNodeGroup;
import etomica.atom.SpeciesAgent;
import etomica.atom.SpeciesRoot;
import etomica.atom.iterator.AtomIteratorSequenceDirectable;
import etomica.atom.iterator.IteratorDirective;
import etomica.junit.UnitTest;


/**
 * Unit test for AtomIteratorSequenceDirectable
 * 
 * @author David Kofke
 *
 */

/*
 * History
 * Created on Jun 6, 2005 by kofke
 */
public class AtomIteratorSequenceDirectableTest extends ListIteratorTest {

    public AtomIteratorSequenceDirectableTest() {
        super();
        UnitTest.VERBOSE = false;
    }

    public void setUp() {
        iterator = new AtomIteratorSequenceDirectable();
    }
    
    /**
     * Iterates over a (no-tab) childlist of a species agent, ensures the
     * appropriate number of iterates are given, when skipping 0 or all, and 
     * when iterating up, down, or both.
     */
    public void testChildListIteration() {
        int nAtoms = 11;
        SpeciesRoot root = UnitTest.makeStandardSpeciesTree(
                null,0,new int[] {nAtoms},null,null);
        SpeciesAgent agent = (SpeciesAgent)((AtomTreeNodeGroup)root.node).getDescendant(new int[] {0,0});
        AtomList list = ((AtomTreeNodeGroup)agent.node).childList;
        int count = list.size();

        //Iterate in different directions and check count of iterates
        //start from first atom
        iterator.setAtom(list.getFirst());
        for(int i=0; i<count+1; i++) {
            iterator.setNumToSkip(i);
            LinkedList[] lists0 = testDirections(iterator);
            testSkip(iterator);
            int upCount = Math.max(0, count-i);
            int dnCount = Math.max(0, 1-i);
            int allCount = upCount + dnCount;
            if(i == 0) allCount--;
            
            assertEquals(lists0[0].size(), upCount);
            assertEquals(lists0[1].size(), dnCount);
            assertEquals(lists0[2].size(), allCount);
        }
        
        //same, starting from last atom
        iterator.setAtom(list.getLast());
        for(int i=0; i<count+1; i++) {
            iterator.setNumToSkip(i);
            LinkedList[] lists0 = testDirections(iterator);
            testSkip(iterator);
            int upCount = Math.max(0, 1-i);
            int dnCount = Math.max(0, count-i);
            int allCount = upCount + dnCount;
            if(i == 0) allCount--;
            
            assertEquals(lists0[0].size(), upCount);
            assertEquals(lists0[1].size(), dnCount);
            assertEquals(lists0[2].size(), allCount);
        }
        
        //same, starting from middle atom
        iterator.setAtom(list.get((count-1)/2));
        for(int i=0; i<count+1; i++) {
            iterator.setNumToSkip(i);
            LinkedList[] lists1 = testDirections(iterator);
            testSkip(iterator);
            int upCount = Math.max(0, (count+1)/2-i);
            int dnCount = Math.max(0, (count+1)/2-i);
            int allCount = upCount + dnCount;
            if(i == 0) allCount--;
            
            assertEquals(lists1[0].size(), upCount);
            assertEquals(lists1[1].size(), dnCount);
            assertEquals(lists1[2].size(), allCount);
        }

        //same, starting from middle atom + 2
        iterator.setAtom(list.get(2+(count-1)/2));
        for(int i=0; i<count+1; i++) {
            iterator.setNumToSkip(i);
            LinkedList[] lists1 = testDirections(iterator);
            testSkip(iterator);
            int upCount = Math.max(0, (count+1)/2-i-2);
            int dnCount = Math.max(0, (count+1)/2-i+2);
            int allCount = upCount + dnCount;
            if(i == 0) allCount--;
            
            assertEquals(lists1[0].size(), upCount);
            assertEquals(lists1[1].size(), dnCount);
            assertEquals(lists1[2].size(), allCount);
        }
    }

    /**
     * Called from testListVariations
     */
    public void iteratorStateTests(AtomList list) {
        iterator.setFirst(list.header.next);
        generalIteratorMethodTests(iterator);

        iterator.setFirst(list.header.previous);
        generalIteratorMethodTests(iterator);
    }
    
    private void testSkip(AtomIteratorSequenceDirectable iterator) {
        for(int i=0; i<11; i++) {
            iterator.setNumToSkip(i);
            testDirections(iterator);
        }
    }
    
    private LinkedList[] testDirections(AtomIteratorSequenceDirectable iterator) {
        iterator.setDirection(IteratorDirective.UP);
        LinkedList listUp = generalIteratorMethodTests(iterator);
        iterator.setDirection(IteratorDirective.DOWN);
        LinkedList listDn = generalIteratorMethodTests(iterator);
        iterator.setDirection(null);
        LinkedList list = generalIteratorMethodTests(iterator);
        return new LinkedList[] {listUp, listDn, list};
    }

    AtomIteratorSequenceDirectable iterator;
}
