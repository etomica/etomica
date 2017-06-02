package etomica.interfacial;

import etomica.potential.PotentialMaster;
import etomica.api.IRandom;
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
        box.removeMolecule(atom.getParentGroup());
        box.addMolecule(atom.getParentGroup());
        return true;
    }
    
    public void rejectNotify() {
        super.rejectNotify();
        // get neighbor manager to notice what we've done
        box.removeMolecule(atom.getParentGroup());
        box.addMolecule(atom.getParentGroup());
    }
}
