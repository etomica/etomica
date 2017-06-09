/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.association;
import etomica.molecule.IMolecule;

/**
 * Interface for a method that defines whether two molecules are considered associated.
 */
 
public interface AssociationDefinitionMolecule {
    
    public boolean isAssociated(IMolecule moleculeA, IMolecule moleculeB);
    
}
