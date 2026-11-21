/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.util.random;

/**
 * Interface for random number generation.
 *
 * @author Andrew Schultz
 */
public interface IRandom {

    /**
     * Convenience method that returns an array of integers, useful as seeds to
     * another RNG.  RandomNumberGeneratorUnix can be instantiated, it is
     * used to return 4 integers.  Otherwise, System.nanotime is split into 2
     * integers.
     */
    static int[] getRandSeedArray() {
        int[] rv;
        if (!RandomNumberGeneratorUnix.hasRandomFile()) {
            // no real problem, probably just not on a unix box.
            rv = new int[2];
            long t = System.nanoTime();
            rv[0] = (int) t; // lower 32 bits
            rv[1] = (int) (t >> 32); // upper 32 bits
        } else {
            RandomNumberGeneratorUnix rng = new RandomNumberGeneratorUnix();
            rv = new int[4];
            for (int i = 0; i < rv.length; i++) {
                rv[i] = rng.nextInt();
            }
            rng.dispose();
        }
        return rv;
    }


    /**
     * Returns a pseudorandom double, uniformly distributed between
     * 0.0 (exclusive) and 1.0 (exclusive).  This is a floating-point value,
     * meaning that the any double can be returned; the double can be
     * arbitrarily small, but will not be 0.
     */
    double nextDouble();

    /**
     * Returns a pseudorandom double, uniformly distributed between
     * 0.0 (inclusive) and 1.0 (exclusive).  For computational speed, this
     * implementation returns a fixed-point random number (at least in binary),
     * meaning that it will be a random integer times 2^-53.
     */
    double nextFixedDouble();

    /**
     * Returns a pseudorandom integer, uniformly distributed between
     * 0 (inclusive) and maxInt (exclusive).
     */
    int nextInt(int maxInt);

    /**
     * Returns a pseudorandom double, Gaussian ("normally") distributed
     * value with mean 0.0 and standard deviation 1.0.
     */
    double nextGaussian();

}
