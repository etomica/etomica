/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom;


/**
 * Interface for a class that screens atoms according
 * to some criterion.
 */
public interface AtomFilter {
    
    /**
     * Returns true if atom is passes test of filter.
     */
    public boolean accept(IAtom a);
    
    public boolean accept(IMolecule mole);

}
