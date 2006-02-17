package etomica.nbr;

import java.util.Iterator;
import java.util.LinkedList;

import etomica.atom.AtomType;
import etomica.atom.iterator.IteratorDirective;
import etomica.nbr.PotentialCalculationUpdateTypeList.PotentialAtomTypeWrapper;
import etomica.phase.Phase;
import etomica.phase.PhaseAgentManager;
import etomica.phase.PhaseAgentManager.PhaseAgentSource;
import etomica.potential.PotentialArray;
import etomica.potential.PotentialMaster;
import etomica.space.Space;

public abstract class PotentialMasterNbr extends PotentialMaster {

    protected PotentialMasterNbr(Space space, PhaseAgentSource phaseAgentSource, 
            PhaseAgentManager phaseAgentManager) {
        super(space);
        this.phaseAgentSource = phaseAgentSource;
        this.phaseAgentManager = phaseAgentManager;
    }

    /**
     * Determines potentials that apply to each AtomType using a special PotentialCalculation
     * which is performed by the superclass (non-neighbor PotentialMaster).
     */
    public void updateTypeList(Phase phase) {
        PotentialCalculationUpdateTypeList pc = new PotentialCalculationUpdateTypeList();
        super.calculate(phase, new IteratorDirective(), pc);
        LinkedList newPotentialTypeList = pc.getPotentialsTypeList();
        Iterator iterator = newPotentialTypeList.iterator();
        while (iterator.hasNext()) {
            addToPotentialTypeList((PotentialAtomTypeWrapper)iterator.next());
        }
    }

    protected void addToPotentialTypeList(PotentialAtomTypeWrapper wrapper) {
        for (int i=0; i<wrapper.atomTypes.length; i++) {
            while (potentialAtomTypeList.length < wrapper.atomTypes[i].getIndex()+1) {
                potentialAtomTypeList = (PotentialArray[])etomica.util.Arrays.addObject(potentialAtomTypeList, new PotentialArray());
            }
            PotentialArray potentialAtomType = potentialAtomTypeList[wrapper.atomTypes[i].getIndex()];
            potentialAtomType.addPotential(wrapper.potential,wrapper.iterator);
            wrapper.atomTypes[i].setInteracting(true);
        }
    }

    public PotentialArray getPotentials(AtomType atomType) {
        return potentialAtomTypeList[atomType.getIndex()];
    }

    public PhaseAgentManager getCellAgentManager() {
        return phaseAgentManager;
    }

    protected PotentialArray[] potentialAtomTypeList = new PotentialArray[0];
    protected PhaseAgentSource phaseAgentSource;
    protected PhaseAgentManager phaseAgentManager;
}
