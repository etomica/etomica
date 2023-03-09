/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.atom;

/**
 * Interface for a class that determines if a pair of atom satisfies a test.
 */
public interface AtomPairTest {

    /**
     * Returns true if the pair of atoms passes the test.
     */
    boolean test(IAtom atom1, IAtom atom2);

}
