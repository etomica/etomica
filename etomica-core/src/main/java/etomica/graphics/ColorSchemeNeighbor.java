/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graphics;

import java.awt.Color;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.simulation.Simulation;
import etomica.nbr.list.NeighborListManager;
import etomica.nbr.list.PotentialMasterList;

/**
 * Color atoms based on being neighbors of the reference atom
 *
 * @author Andrew Schultz
 */
public class ColorSchemeNeighbor extends ColorSchemeCollectiveAgent {
    
    public ColorSchemeNeighbor(Simulation sim, PotentialMasterList potentialMaster, Box box) {
        super(box);
        typeColorScheme = new ColorSchemeByType(sim);
        leafList = box.getLeafList();
        neighborManager = potentialMaster.getNeighborManager(box);
    }
    
    public void colorAllAtoms() {
		//color all atoms according to their type
        int nLeaf = leafList.getAtomCount();
        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            IAtom atom = leafList.getAtom(iLeaf);
            agentManager.setAgent(atom, typeColorScheme.getAtomColor(atom));
        }
        if (referenceAtom == null) {
            return;
        }
        //color blue the neighbor atoms in same group
        IAtomList[] list = neighborManager.getDownList(referenceAtom);
        for (int foo=0; foo<2; foo++) {
            if (foo==1) {
                list = neighborManager.getUpList(referenceAtom);
            }
            for (int i=0; i<list.length; i++) {
                for (int j=0; j<list[i].getAtomCount(); j++) {
                    IAtom atom = list[i].getAtom(j);
                    if(atom.getType() == referenceAtom.getType()) {
                        agentManager.setAgent(atom, Color.blue);
                    } else {
                        agentManager.setAgent(atom, Color.yellow);
                    }
                }
            }
        }
        //color green the target atom 
        agentManager.setAgent(referenceAtom, Color.green);
    }
    
    public void setAtom(IAtom a) {
        referenceAtom = a;
    }

    public IAtom getAtom() {
        return referenceAtom;
    }
    
    private static final long serialVersionUID = 1L;
    private IAtom referenceAtom;
    private final NeighborListManager neighborManager;
    private final IAtomList leafList;
    private final ColorSchemeByType typeColorScheme;
}
