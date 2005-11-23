package etomica.nbr;

import java.util.Iterator;
import java.util.LinkedList;

import etomica.atom.AtomType;
import etomica.atom.iterator.IteratorDirective;
import etomica.atom.iterator.IteratorFactory;
import etomica.nbr.PotentialCalculationUpdateTypeList.PotentialAtomTypeWrapper;
import etomica.phase.Phase;
import etomica.potential.Potential;
import etomica.potential.PotentialArray;
import etomica.potential.PotentialMaster;
import etomica.space.Space;

public class PotentialMasterNbr extends PotentialMaster {

    protected PotentialMasterNbr(Space space) {
        super(space);
    }

    protected PotentialMasterNbr(Space space, IteratorFactory iteratorFactory) {
        super(space, iteratorFactory);
    }

    /**
     * Performs cell-assignment potentialCalculation.  Assigns all molecules
     * to their cells, and invokes superclass method causing setup to be
     * performed iterating using species/potential hierarchy.
     */
    public void updateTypeList(Phase phase) {
        PotentialCalculationUpdateTypeList pc = new PotentialCalculationUpdateTypeList();
        super.calculate(phase, new IteratorDirective(), pc);
        LinkedList newPotentialTypeList = pc.getPotentialsTypeList();
        Iterator iterator = newPotentialTypeList.iterator();
        while (iterator.hasNext()) {
            PotentialAtomTypeWrapper wrapper = (PotentialAtomTypeWrapper)iterator.next();
            addToPotentialTypeList(wrapper.potential,wrapper.atomTypes);
        }
    }

    protected void addToPotentialTypeList(Potential potential, AtomType[] atomTypes) {
        for (int i=0; i<atomTypes.length; i++) {
            while (potentialAtomTypeList.length < atomTypes[i].getIndex()+1) {
                potentialAtomTypeList = (PotentialArray[])etomica.util.Arrays.addObject(potentialAtomTypeList, new PotentialArray());
            }
            PotentialArray potentialAtomType = potentialAtomTypeList[atomTypes[i].getIndex()];
            potentialAtomType.addPotential(potential);
            atomTypes[i].setInteracting(true);
        }
    }

    public PotentialArray getPotentials(AtomType atomType) {
        return potentialAtomTypeList[atomType.getIndex()];
    }

    protected PotentialArray[] potentialAtomTypeList = new PotentialArray[0];
}
