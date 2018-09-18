/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.math.discrete;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;

/**
 * Created by alex on 5/1/17.
 */
public class FixedSumIteratorTest {
    private FixedSumIterator iter;

    @BeforeEach
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
