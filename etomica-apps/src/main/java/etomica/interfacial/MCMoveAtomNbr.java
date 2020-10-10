/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.interfacial;

import etomica.potential.PotentialMaster;
import etomica.space.Vector;
import etomica.util.random.IRandom;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.space.Space;

public class MCMoveAtomNbr extends MCMoveAtom {

    public MCMoveAtomNbr(IRandom random, PotentialMaster potentialMaster, Space _space) {
        super(random, potentialMaster, _space);
    }

    public boolean doTrial() {
        boolean rv = super.doTrial();
        if (!rv) return false;
        // get neighbor manager to notice what we've done
        Vector pos = atom.getPosition().duplicate();
        box.removeMolecule(atom.getParentGroup());
        box.addNewMolecule(atom.getParentGroup().getType(), mol -> mol.getChildList().get(0).getPosition().E(pos));
        return true;
    }
    
    public void rejectNotify() {
        super.rejectNotify();
        // get neighbor manager to notice what we've done
        Vector pos = atom.getPosition().duplicate();
        box.removeMolecule(atom.getParentGroup());
        box.addNewMolecule(atom.getParentGroup().getType(), mol -> mol.getChildList().get(0).getPosition().E(pos));
    }
}
