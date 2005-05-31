package etomica.junit.atom.iterator;

import etomica.atom.iterator.ApiListSimple;
import etomica.atom.iterator.AtomsetIteratorListDependent;


/**
 * Performs unit tests of the ApiListSimple iterator
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
    }
    
    /**
     * Tests iterator as given (from testListVariations method) using just
     * generalIteratorMethodTests.  
     * Test of correctness of iterates is concluded from agreement between hasNext/next
     * and allAtoms iterations, which is among the tests done here.
     */
    public void iteratorStateTests(AtomsetIteratorListDependent iterator) {
        generalIteratorMethodTests(iterator);
    }
}
