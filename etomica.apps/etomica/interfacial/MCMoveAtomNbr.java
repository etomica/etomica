package etomica.interfacial;

import etomica.api.IPotentialMaster;
import etomica.api.IRandom;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.space.ISpace;

public class MCMoveAtomNbr extends MCMoveAtom {

    public MCMoveAtomNbr(IRandom random, IPotentialMaster potentialMaster, ISpace _space) {
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
