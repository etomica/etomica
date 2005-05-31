/*
 * Created on Oct 7, 2004
 */
package etomica.junit.atom.iterator;

import etomica.atom.iterator.AtomsetIteratorListDependent;
import etomica.atom.iterator.AtomIteratorListSimple;
import etomica.junit.UnitTest;



/**
 * Class for testing simple lists using JUNIT.
 * 
 * @author Ken Benjamin
 */
public class ListIteratorTestSimple extends ListIteratorTest {

	public ListIteratorTestSimple() {
		super(new AtomIteratorListSimple());
        UnitTest.VERBOSE = false;
	}
/**
 * does nothing 
 */
	public void setUp() {
		//does nothing
	}
	
/**
 * Calls the generalIteratorMethodTests method, which takes the iterator and
 * performs tests on all the general methods contained in the class
 * AtomIteratorListSimple.
 */	
	public void iteratorStateTests(AtomsetIteratorListDependent iterator) {
		generalIteratorMethodTests(iterator);
	}
		

}



