package etomica.nbr;

import etomica.atom.AtomType;
import etomica.phase.PhaseAgentManager;
import etomica.phase.PhaseAgentManager.PhaseAgentSource;
import etomica.potential.Potential;
import etomica.potential.PotentialArray;
import etomica.potential.PotentialGroup;
import etomica.potential.PotentialMaster;
import etomica.space.Space;
import etomica.species.Species;
import etomica.util.Arrays;

public abstract class PotentialMasterNbr extends PotentialMaster {

    protected PotentialMasterNbr(Space space, PhaseAgentSource phaseAgentSource, 
            PhaseAgentManager phaseAgentManager) {
        super(space);
        this.phaseAgentSource = phaseAgentSource;
        this.phaseAgentManager = phaseAgentManager;
    }
    
    public PotentialGroup makePotentialGroup(int nBody) {
        return new PotentialGroupNbr(nBody, space);
    }
    
    public void addPotential(Potential potential, Species[] species) {
        super.addPotential(potential, species);
        if (!(potential instanceof PotentialGroup)) {
             AtomType[] atomTypes = moleculeTypes(species);
             if (potential.getRange() == Double.POSITIVE_INFINITY) {
                 System.err.println("You gave me a molecular range-independent potential and I'm very confused now");
                 return;
             }
             //the potential is range-dependent 
             for (int i=0; i<atomTypes.length; i++) {
                 addRangedPotential(potential,atomTypes[i]);
             }
        }
    }
    
    public void potentialAddedNotify(Potential subPotential, PotentialGroup pGroup) {
        super.potentialAddedNotify(subPotential, pGroup);
        AtomType[] atomTypes = pGroup.getAtomTypes(subPotential);
        if (atomTypes == null) {
            if (pGroup.nBody() == 1 && subPotential.getRange() == Double.POSITIVE_INFINITY) {
                boolean found = false;
                for (int i=0; i<allPotentials.length; i++) {
                    if (allPotentials[i] == pGroup) {
                        found = true;
                    }
                }
                if (!found) {
                    allPotentials = (Potential[])etomica.util.Arrays.addObject(allPotentials, pGroup);
                }
                //pGroup is PotentialGroupNbr
                AtomType[] parentType = getAtomTypes(pGroup);
                while (parentType[0].getIndex() > rangedPotentialAtomTypeList.length-1) {
                    rangedPotentialAtomTypeList = (PotentialArray[])etomica.util.Arrays.addObject(rangedPotentialAtomTypeList, new PotentialArray());
                    intraPotentialAtomTypeList = (PotentialArray[])etomica.util.Arrays.addObject(intraPotentialAtomTypeList, new PotentialArray());
                }
                intraPotentialAtomTypeList[parentType[0].getIndex()].addPotential(pGroup,null);
            }
            else {
                //FIXME what to do with this case?  Fail!
                System.err.println("You have a child-potential of a 2-body PotentialGroup or range-dependent potential, but it's not type-based.  Enjoy crashing or fix bug 85");
            }
            return;
        }
        if (subPotential.getRange() == Double.POSITIVE_INFINITY) {
            //FIXME what to do with this case?
            System.err.println("you have an infinite-ranged potential that's type based!  I don't like you.");
            return;
        }
        for (int i=0; i<atomTypes.length; i++) {
            addRangedPotential(subPotential,atomTypes[i]);
        }
    }
    
    protected void addRangedPotential(Potential potential, AtomType atomType) {
        while (rangedPotentialAtomTypeList.length < atomType.getIndex()+1) {
            rangedPotentialAtomTypeList = (PotentialArray[])etomica.util.Arrays.addObject(rangedPotentialAtomTypeList, new PotentialArray());
            intraPotentialAtomTypeList = (PotentialArray[])etomica.util.Arrays.addObject(intraPotentialAtomTypeList, new PotentialArray());
        }
        PotentialArray potentialAtomType = rangedPotentialAtomTypeList[atomType.getIndex()];
        potentialAtomType.addPotential(potential,null);
        atomType.setInteracting(true);
        boolean found = false;
        for (int i=0; i<allPotentials.length; i++) {
            if (allPotentials[i] == potential) {
                found = true;
            }
        }
        if (!found) {
            allPotentials = (Potential[])etomica.util.Arrays.addObject(allPotentials, potential);
        }
    }
    
    public void removePotential(Potential potential) {
        super.removePotential(potential);
        if (potential.getRange() < Double.POSITIVE_INFINITY) {
            rangedPotentialAtomTypeList = (PotentialArray[])Arrays.removeObject(rangedPotentialAtomTypeList,potential);
        }
        else if (potential instanceof PotentialGroup) {
            intraPotentialAtomTypeList = (PotentialArray[])Arrays.removeObject(intraPotentialAtomTypeList,potential);
        }
        allPotentials = (Potential[])Arrays.removeObject(allPotentials,potential);
    }
    
    public PotentialArray getRangedPotentials(AtomType atomType) {
        return rangedPotentialAtomTypeList[atomType.getIndex()];
    }

    public PotentialArray getIntraPotentials(AtomType atomType) {
        return intraPotentialAtomTypeList[atomType.getIndex()];
    }
    
    public PhaseAgentManager getCellAgentManager() {
        return phaseAgentManager;
    }

    protected PotentialArray[] rangedPotentialAtomTypeList = new PotentialArray[0];
    protected PotentialArray[] intraPotentialAtomTypeList = new PotentialArray[0];
    protected Potential[] allPotentials = new Potential[0];
    protected PhaseAgentSource phaseAgentSource;
    protected PhaseAgentManager phaseAgentManager;
}
