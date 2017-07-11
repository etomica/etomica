/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.association;

import etomica.atom.IAtom;
import etomica.atom.AtomArrayList;

/**
 * Interface for a class that is capable of populating a list of atoms in an
 * smer as well as checking the validity of the bonding 
 * @author andrew
 *
 */
public interface IAssociationHelper {

    /**
     * Populates smerList with all atoms in the smer that contains atom.
     * If mightBeBroken is false, an exception is throw if any invalid bonding
     * is encountered (along with potentially useful information about the
     * invalid bonding).  If mightBeBroken is true, then populateList returns 
     * true if any invalid bonding is encountered.  populateList returns false
     * if invalid bonding is not encountered.
     */
    public boolean populateList(AtomArrayList smerList, IAtom atom,
            boolean mightBeBroken);

}
