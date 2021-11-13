/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.potential;

/**
 * Factory for a truncated potential.
 */
public interface TruncationFactory {

    /**
     * @return a truncated potential that wraps (and sums) the given pair
     * potentials.
     */
    Potential2Soft make(Potential2Soft... p2);
}
