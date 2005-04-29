/*
 * Created on Oct 7, 2004
 *
 * TODO To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package etomica.junit;
//import java.util.Collections;
//import java.util.LinkedList;

import etomica.*;
import etomica.atom.iterator.AtomIteratorListDependent;
import etomica.atom.iterator.AtomIteratorListSimple;
import junit.framework.*;
import java.util.*;



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
	public void iteratorStateTests(AtomIteratorListDependent iterator) {
		generalIteratorMethodTests(iterator);
	}
		

}



