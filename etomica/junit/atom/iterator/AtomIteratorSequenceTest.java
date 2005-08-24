package etomica.junit.atom.iterator;

import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedList;

import etomica.IteratorDirective;
import etomica.action.AtomsetAction;
import etomica.atom.Atom;
import etomica.atom.AtomLinker;
import etomica.atom.AtomList;
import etomica.atom.AtomTreeNodeGroup;
import etomica.atom.SpeciesAgent;
import etomica.atom.SpeciesRoot;
import etomica.atom.iterator.AtomIteratorSequence;
import etomica.junit.UnitTest;

/**
 * Unit test for AtomIteratorSequencer
 */

/*
 * Created on June 6, 2005 by kofke
 */
public class AtomIteratorSequenceTest extends ListIteratorTest {

    public AtomIteratorSequenceTest() {
        super();
        UnitTest.VERBOSE = false;
    }

    public void setUp() {
        upIterator = new AtomIteratorSequence(IteratorDirective.UP);
        dnIterator = new AtomIteratorSequence(IteratorDirective.DOWN);

    }
    
    /**
     * Iterates over a (no-tab) childlist of a species agent, ensures the
     * appropriate number of iterates are given, and that up and down give
     * list in opposite orders.
     */
    public void testChildListIteration() {
        int nAtoms = 10;
        SpeciesRoot root = UnitTest.makeStandardSpeciesTree(
                null,0,new int[] {nAtoms},null,null);
        SpeciesAgent agent = (SpeciesAgent)((AtomTreeNodeGroup)root.node).getDescendant(new int[] {0,0});
        AtomList list = ((AtomTreeNodeGroup)agent.node).childList;
        upIterator.setAtom(list.getFirst());
        dnIterator.setAtom(list.getFirst());
        LinkedList listUp0 = generalIteratorMethodTests(upIterator);
        LinkedList listDn0 = generalIteratorMethodTests(dnIterator);
        upIterator.setAtom(list.getLast());
        dnIterator.setAtom(list.getLast());
        LinkedList listUp1 = generalIteratorMethodTests(upIterator);
        LinkedList listDn1 = generalIteratorMethodTests(dnIterator);
        assertEquals(listUp0.size(),nAtoms);
        assertEquals(listDn0.size(), 1);
        assertEquals(listUp1.size(), 1);
        assertEquals(listDn1.size(),nAtoms);
        Collections.reverse(listDn1);
        Object[] list0 = listUp0.toArray();
        Object[] list1 = listDn1.toArray();
        assertTrue(Arrays.equals(list0, list1));
        
        //check that allAtoms doesn't clobber iteration state
        //(as documented for this iterator)
        upIterator.setAtom(list.getFirst());
        upIterator.reset();
        for(int i=0; i<nAtoms/2; i++) upIterator.nextAtom();
        Atom next = (Atom)upIterator.peek();
        upIterator.allAtoms(AtomsetAction.NULL);
        assertTrue(upIterator.hasNext());
        assertEquals(next, upIterator.next());
    }
    
    /**
     * Performs generalIteratorMethodTests method for the given list,
     * for up-list and down-list iterators, beginning at the entry
     * after the header, and the last entry before the header.  Performs
     * no tests to examine the specific iterates.
     */
    public void iteratorStateTests(AtomList list) {
        upIterator.setFirst(list.header.next);
        dnIterator.setFirst(list.header.next);
        generalIteratorMethodTests(upIterator);
        generalIteratorMethodTests(dnIterator);

        checkNextLinker(upIterator, dnIterator);

        upIterator.setFirst(list.header.previous);
        dnIterator.setFirst(list.header.previous);
        generalIteratorMethodTests(upIterator);
        generalIteratorMethodTests(dnIterator);
        
        checkNextLinker(upIterator, dnIterator);
    }
    
    private void checkNextLinker(AtomIteratorSequence upIterator, AtomIteratorSequence dnIterator) {
        AtomLinker first = upIterator.getFirst();
        assertEquals(first, dnIterator.getFirst());
        upIterator.reset();
        int count = 0;
        //loop through upIterator until return to first
        while(count < 1000000) {//ensure no infinite loop
            if((upIterator.nextLinker() == first) && (count != 0)) break;//check of count must be 2nd part of &&
            count++;
        }
        print("Count: "+count);
        //loop through dnIterator same number of times and see that first is recovered
        dnIterator.reset();
        for(int i=0; i<count; i++) {
            dnIterator.nextLinker();
        }
        assertEquals(dnIterator.nextLinker(), first);
    }

    private AtomIteratorSequence upIterator, dnIterator;
}

