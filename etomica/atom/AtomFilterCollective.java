/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom;

/**
 * Filter that needs to be informed whenever something has changed and the
 * filtered status of atoms should be redetermined.
 * 
 * @author Andrew Schultz
 */
public interface AtomFilterCollective extends AtomFilter {

    /**
     * Informs the filter that something in the system has changed (atom
     * positions, velocityies, etc) has changed and the filtered status of the
     * atoms should be determined.
     */
    public void resetFilter();
}
