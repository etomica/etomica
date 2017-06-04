/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graphics;

import java.awt.Color;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.simulation.Simulation;
import etomica.nbr.PotentialMasterNbr;
import etomica.nbr.cell.Api1ACell;

/**
 * Colors atoms based on being in a neighboring cell to a reference atom.
 */
public class ColorSchemeNeighborCell extends ColorSchemeCollectiveAgent {
    
    public ColorSchemeNeighborCell(Simulation sim, PotentialMasterNbr potentialMaster, Box box, int dim) {
        super(box);
        typeColorScheme = new ColorSchemeByType(sim);
        leafList = box.getLeafList();
        nbrIterator = new Api1ACell(dim, 1.0, potentialMaster.getCellAgentManager());
        nbrIterator.setDirection(null);
        nbrIterator.setBox(box);
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
        nbrIterator.reset();
        for (IAtomList pair = nbrIterator.next(); pair != null;
             pair = nbrIterator.next()) {
            IAtom atom = pair.getAtom(1);
            if(atom.getType() == referenceAtom.getType()) {
                agentManager.setAgent(atom, Color.blue);
            } else {
                agentManager.setAgent(atom, Color.yellow);
            }
        }
        //color green the target atom 
        agentManager.setAgent(referenceAtom, Color.green);
    }
    
    public void setAtom(IAtom a) {
        referenceAtom = a;
        nbrIterator.setTarget(a);
    }

    public IAtom getAtom() {
        return referenceAtom;
    }
    
    private static final long serialVersionUID = 1L;
    private IAtom referenceAtom;
    private final Api1ACell nbrIterator;
    private final IAtomList leafList;
    private final ColorSchemeByType typeColorScheme;
}
