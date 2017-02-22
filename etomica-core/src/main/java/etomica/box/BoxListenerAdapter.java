/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.box;

import etomica.api.IBoxAtomEvent;
import etomica.api.IBoxAtomIndexEvent;
import etomica.api.IBoxEvent;
import etomica.api.IBoxIndexEvent;
import etomica.api.IBoxListener;
import etomica.api.IBoxMoleculeCountEvent;
import etomica.api.IBoxMoleculeEvent;
import etomica.api.IBoxMoleculeIndexEvent;

public class BoxListenerAdapter implements IBoxListener {

    public void boxAtomAdded(IBoxAtomEvent e) { }
    
    public void boxAtomRemoved(IBoxAtomEvent e) { }
    
    public void boxMoleculeAdded(IBoxMoleculeEvent e) { }
    
    public void boxMoleculeRemoved(IBoxMoleculeEvent e) { }
    
    public void boxGlobalAtomIndexChanged(IBoxIndexEvent e) { }
    
    public void boxGlobalAtomLeafIndexChanged(IBoxIndexEvent e) { }
    
    public void boxAtomLeafIndexChanged(IBoxAtomIndexEvent e) { }
    
    public void boxMoleculeIndexChanged(IBoxMoleculeIndexEvent e) { }
    
    public void boxNumberMolecules(IBoxMoleculeCountEvent e) { }
    
    public void boxNeighborsUpdated(IBoxEvent e) { }
    
}
