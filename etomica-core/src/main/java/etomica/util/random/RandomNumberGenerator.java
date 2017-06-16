/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.util.random;

import java.util.Random;

/**
 * Random number generator class.  This just wraps Java's own RNG.
 *
 * @author Andrew Schultz
 */
public class RandomNumberGenerator implements IRandom {

    private final Random jRandom;

    public RandomNumberGenerator() {
        jRandom = new Random();
    }

    public RandomNumberGenerator(long seed) {
        jRandom = new Random(seed);
    }

    // XXX we're just wrapping Java's implementation, which will be fixed-point
    public double nextDouble() {
        return jRandom.nextDouble();
    }

    public double nextFixedDouble() {
        return jRandom.nextDouble();
    }

    public int nextInt(int maxInt) {
        return jRandom.nextInt(maxInt);
    }

    public double nextGaussian() {
        return jRandom.nextGaussian();
    }

    public Random getWrappedRandom() {
        return jRandom;
    }
}
