/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.api;

public interface IBoxListener {

    public void boxAtomAdded(IBoxAtomEvent e);
    
    public void boxAtomRemoved(IBoxAtomEvent e);
    
    public void boxMoleculeAdded(IBoxMoleculeEvent e);
    
    public void boxMoleculeRemoved(IBoxMoleculeEvent e);
    
    public void boxGlobalAtomLeafIndexChanged(IBoxIndexEvent e);
    
    public void boxGlobalAtomIndexChanged(IBoxIndexEvent e);
    
    public void boxAtomLeafIndexChanged(IBoxAtomIndexEvent e);
    
    public void boxMoleculeIndexChanged(IBoxMoleculeIndexEvent e);
    
    public void boxNumberMolecules(IBoxMoleculeCountEvent e);

}
