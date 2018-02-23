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

public class ColorSchemeCoordination extends ColorSchemeCollectiveAgent {
    
    public ColorSchemeCoordination(Simulation sim, PotentialMasterList potentialMaster, Box box) {
        super(box);

        leafList = box.getLeafList();
        nbrManager = potentialMaster.getNeighborManager(box);
        colors = new Color[]{Color.WHITE,Color.MAGENTA, Color.RED, Color.ORANGE,Color.GRAY,Color.GREEN,Color.BLUE,Color.WHITE};
    }
    
    public void colorAllAtoms() {
    	
		//color all atoms according to their type
        int nLeaf = leafList.size();
        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            IAtom atom = leafList.get(iLeaf);
            
            int coordNum = nbrManager.getDownList(atom)[0].size();
            coordNum += nbrManager.getUpList(atom)[0].size();
            
            if(coordNum>colors.length-1){
            	coordNum = colors.length-1;	
            }
            agentManager.setAgent(atom, colors[coordNum]);
        }
    }
        
    private static final long serialVersionUID = 1L;
    private final IAtomList leafList;
    private Color[] colors;
    private NeighborListManager nbrManager;
}
