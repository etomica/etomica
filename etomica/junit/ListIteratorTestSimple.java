/*
 * Created on Oct 7, 2004
 */
package etomica.junit;

import etomica.atom.iterator.AtomsetIteratorListDependent;
import etomica.atom.iterator.AtomIteratorListSimple;



/**
 * @author Ken Benjamin
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 * 
 * Class for testing simple lists using JUNIT.
 */
public class ListIteratorTestSimple extends ListIteratorTest {

	public ListIteratorTestSimple() {
		super(new AtomIteratorListSimple());
        UnitTest.VERBOSE = true;
	}
/**
 * setUp is a required method for any JUNIT test.  Here, a new iterator of type
 * AtomIteratorListSimple is constructed for subsequent use in unit testing.		 
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
		IteratorTest.generalIteratorMethodTests(iterator);
	}
		

}



