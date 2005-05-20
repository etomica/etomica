package etomica.junit;

import etomica.atom.iterator.ApiListSimple;
import etomica.atom.iterator.AtomsetIteratorListDependent;


/**
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 *
 * @author David Kofke
 *
 */

/*
 * History
 * Created on May 19, 2005 by kofke
 */
public class ApiListSimpleTest extends ListIteratorTest {

    /**
     * 
     */
    public ApiListSimpleTest() {
        super(new ApiListSimple());
        UnitTest.VERBOSE = true;
    }
    

    public void iteratorStateTests(AtomsetIteratorListDependent iterator) {
        IteratorTest.generalIteratorMethodTests(iterator);
    }
}
