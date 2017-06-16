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
