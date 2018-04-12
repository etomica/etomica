/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graphics;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.nbr.PotentialMasterNbr;
import etomica.nbr.cell.Api1ACell;
import etomica.nbr.cell.NeighborCellManager;

import java.awt.*;

/**
 * Colors atoms based on being in a neighboring cell to a reference atom.
 */
public class ColorSchemeNeighborCell extends ColorSchemeCollectiveAgent {

    private final Api1ACell nbrIterator;
    private final IAtomList leafList;
    private final ColorSchemeByType typeColorScheme;
    private IAtom referenceAtom;

    public ColorSchemeNeighborCell(PotentialMasterNbr potentialMaster, Box box, int dim) {
        super(box);
        typeColorScheme = new ColorSchemeByType();
        leafList = box.getLeafList();
        nbrIterator = new Api1ACell(1.0, box, (NeighborCellManager) potentialMaster.getBoxCellManager(box));
        nbrIterator.setDirection(null);
    }

    public void colorAllAtoms() {
        //color all atoms according to their type
        int nLeaf = leafList.size();
        for (int iLeaf = 0; iLeaf < nLeaf; iLeaf++) {
            IAtom atom = leafList.get(iLeaf);
            agentManager.setAgent(atom, typeColorScheme.getAtomColor(atom));
        }
        if (referenceAtom == null) {
            return;
        }
        //color blue the neighbor atoms in same group
        nbrIterator.reset();
        for (IAtomList pair = nbrIterator.next(); pair != null;
             pair = nbrIterator.next()) {
            IAtom atom = pair.get(1);
            if (atom.getType() == referenceAtom.getType()) {
                agentManager.setAgent(atom, Color.blue);
            } else {
                agentManager.setAgent(atom, Color.yellow);
            }
        }
        //color green the target atom
        agentManager.setAgent(referenceAtom, Color.green);
    }

    public IAtom getAtom() {
        return referenceAtom;
    }

    public void setAtom(IAtom a) {
        referenceAtom = a;
        nbrIterator.setTarget(a);
    }
}
