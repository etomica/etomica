/*
 * Created on Oct 7, 2004
 */
package etomica.junit.atom.iterator;

import etomica.atom.AtomList;
import etomica.atom.iterator.AtomIteratorListSimple;
import etomica.junit.UnitTest;

/**
 * Unit test for AtomIteratorListSimple.
 * 
 * @author Ken Benjamin
 */
public class ListIteratorTestSimple extends ListIteratorTest {

    public ListIteratorTestSimple() {
        super();
        UnitTest.VERBOSE = false;
    }

    /**
     * Does nothing
     */
    public void setUp() {
        //does nothing
    }

    /**
     * Calls the generalIteratorMethodTests method, which takes the iterator and
     * performs tests on all the general methods contained in the class
     * AtomIteratorListSimple.
     */
    public void iteratorStateTests(AtomList list) {
        generalIteratorMethodTests(new AtomIteratorListSimple(list));
    }

}

