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
        AtomType[] atomTypes = moleculeTypes(species);
        if (potential.getCriterion().isRangeDependent()) {
            //the potential might still be range-independent (intermolecular potential group) 
            for (int i=0; i<atomTypes.length; i++) {
                addRangedPotential(potential,atomTypes[i]);
            }
        }
        else if (!(potential instanceof PotentialGroup)) {
            System.err.println("You gave me a concrete molecular potential range-independent potential and I'm very confused now");
        }
    }
    
    public void potentialAddedNotify(Potential subPotential, PotentialGroup pGroup) {
        super.potentialAddedNotify(subPotential, pGroup);
        AtomType[] atomTypes = pGroup.getAtomTypes(subPotential);
        if (atomTypes == null) {
            if (pGroup.nBody() == 1 && !subPotential.getCriterion().isRangeDependent()) {
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
                PotentialArray potentialAtomType = intraPotentialAtomTypeList[parentType[0].getIndex()];
                potentialAtomType.addPotential(pGroup,null);
            }
            // potential is not type-based
            return;
        }
        if (!subPotential.getCriterion().isRangeDependent()) {
            //FIXME what to do with this case?
            System.err.println("you have a range-dependent potential that's not type based!  I don't like you");
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
