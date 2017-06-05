/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.association;
import etomica.atom.IAtom;
import etomica.space.Space;

public abstract class BiasVolume implements AssociationDefinition, java.io.Serializable {
    
    public final Space space;
    public BiasVolume(Space s) {
        space = s;
    }
    
    public abstract double biasVolume();  
    public abstract void biasInsert(IAtom atomA, IAtom atomB);
    public abstract boolean isAssociated(IAtom atomA, IAtom atomB);
 
    
    
}
