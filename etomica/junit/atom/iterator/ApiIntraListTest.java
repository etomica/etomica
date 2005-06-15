package etomica.junit.atom.iterator;

import etomica.atom.AtomList;
import etomica.atom.iterator.ApiIntraList;


/**
 * Performs unit tests of the ApiIntraList iterator
 *
 * @author David Kofke
 *
 */

/*
 * History
 * Created on May 19, 2005 by kofke
 */
public class ApiIntraListTest extends ListIteratorTest {

    /**
     * 
     */
    public ApiIntraListTest() {
        super();
    }
    
    /**
     * Tests iterator for given list (from testListVariations method) using just
     * generalIteratorMethodTests.  
     * Test of correctness of iterates is concluded from agreement between hasNext/next
     * and allAtoms iterations, which is among the tests done here.
     */
    public void iteratorStateTests(AtomList list) {
        generalIteratorMethodTests(new ApiIntraList(list));
    }
}
