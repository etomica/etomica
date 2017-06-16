/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.association;
import etomica.atom.IAtom;

/**
 * Interface for a method that defines whether two atoms are considered associated.
 */
 
public interface AssociationDefinition {
    
    public boolean isAssociated(IAtom atomA, IAtom atomB);
    
}
