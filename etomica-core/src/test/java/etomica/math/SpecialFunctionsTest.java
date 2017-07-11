package etomica.math;

import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.ExpectedException;

import static org.junit.Assert.*;

/**
 * Created by alex on 4/24/17.
 */
public class SpecialFunctionsTest {
    private static final double EPSILON = 1e-10;

    @Rule
    public ExpectedException thrown = ExpectedException.none();

    @Test
    public void testFactorial() {
        assertEquals(1, SpecialFunctions.factorial(0));
        assertEquals(1, SpecialFunctions.factorial(1));
        assertEquals(120, SpecialFunctions.factorial(5));
        assertEquals(2432902008176640000L, SpecialFunctions.factorial(20));

    }

    @Test
    public void testFactorialException() {
        // this asserts that an exception will be thrown for -1!
        thrown.expect(IllegalArgumentException.class);
        SpecialFunctions.factorial(-1);
    }

    @Test
    public void testLnFactorial() {
        assertEquals(0, SpecialFunctions.lnFactorial(1), EPSILON);
        assertEquals(10.60460290274525, SpecialFunctions.lnFactorial(8), EPSILON);
        assertEquals(42.335616460753485, SpecialFunctions.lnFactorial(20), EPSILON);
        assertEquals(363.73937555556349, SpecialFunctions.lnFactorial(100), EPSILON);
        assertEquals(82108.92783681436, SpecialFunctions.lnFactorial(10000),1e-9);
    }


}
