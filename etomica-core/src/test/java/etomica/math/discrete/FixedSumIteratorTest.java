package etomica.math.discrete;

import org.junit.Before;
import org.junit.Test;

import static org.junit.Assert.*;

/**
 * Created by alex on 5/1/17.
 */
public class FixedSumIteratorTest {
    private FixedSumIterator iter;

    @Before
    public void setUp() throws Exception {
        iter = new FixedSumIterator(5);
        iter.setSum(5);
        iter.reset();
    }

    @Test
    public void testIteratorResults() {
        assertArrayEquals(new int[]{0, 0, 0, 0, 1}, iter.next());
        assertArrayEquals(new int[]{1, 0, 0, 1, 0}, iter.next());
        assertArrayEquals(new int[]{0, 1, 1, 0, 0}, iter.next());
        assertArrayEquals(new int[]{2, 0, 1, 0, 0}, iter.next());
        assertArrayEquals(new int[]{1, 2, 0, 0, 0}, iter.next());
        assertArrayEquals(new int[]{3, 1, 0, 0, 0}, iter.next());
        assertArrayEquals(new int[]{5, 0, 0, 0, 0}, iter.next());
    }

}
