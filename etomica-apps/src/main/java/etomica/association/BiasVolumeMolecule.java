/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.association;

import etomica.atom.IMolecule;
import etomica.space.Space;

public abstract class BiasVolumeMolecule implements AssociationDefinitionMolecule, java.io.Serializable {
    
    public final Space space;
    public BiasVolumeMolecule(Space s) {
        space = s;
    }
    
    public abstract double biasVolume();  
    public abstract void biasInsert(IMolecule moleculeA, IMolecule moleculeB);
    public abstract boolean isAssociated(IMolecule moleculeA, IMolecule moleculeB);
 
    
    
}
