/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.molecule;

import java.io.Serializable;

/**
 * @author Tai Boon Tan  
 * 
 * 
 */
public class MoleculeToMoleculeSetFixed implements MoleculeToMoleculeList, MoleculeToIndex, Serializable {

    private static final long serialVersionUID = 1L;

    public MoleculeToMoleculeSetFixed() {
        moleculeArrayList = new MoleculeArrayList();
    }
    
    public void setArrayList(MoleculeArrayList list) {
        moleculeArrayList = list;
    }
    
    public IMoleculeList getMoleculeList(IMolecule molecule) {
        return moleculeArrayList;
    }
    
    public int getIndex(IMolecule molecule) {
        return moleculeArrayList.indexOf(molecule);
    }

    private MoleculeArrayList moleculeArrayList;
}
