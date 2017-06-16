/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graphics;

import etomica.atom.IAtomList;

/**
 * The BondManager is responsible for creation and disposal of bonds between
 * atom pairs.
 *
 * @author Andrew Schultz
 */
public interface BondManager {
    
    /**
     * Creates a bond object for the given pair with type determined by the
     * bondType object.  The returned object can be used later to
     * tell the BondManager to remove the bond.
     */
    public Object makeBond(IAtomList pair, Object bondType);
    
    /**
     * Notifies the BondManager that the given bond no longer exists in the
     * system.  The given bond object must have been returned from a previous
     * call to makeBond.
     */
    public void releaseBond(Object bond);
}
