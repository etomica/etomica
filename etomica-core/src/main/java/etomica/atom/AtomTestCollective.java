/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom;

/**
 * AtomTest that needs to be informed whenever something has changed and the
 * status of atoms should be redetermined.
 * 
 * @author Andrew Schultz
 */
public interface AtomTestCollective extends AtomTest {

    /**
     * Informs the AtomTest that something in the system has changed (atom
     * positions, velocityies, etc) has changed and the status of the
     * atoms should be determined.
     */
    public void resetTest();
}
