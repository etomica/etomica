package etomica.nbr;

import etomica.atom.AtomType;
import etomica.atom.AtomTypeAgentManager;
import etomica.phase.PhaseAgentManager;
import etomica.phase.PhaseAgentManager.PhaseAgentSource;
import etomica.potential.Potential;
import etomica.potential.PotentialArray;
import etomica.potential.PotentialGroup;
import etomica.potential.PotentialMaster;
import etomica.space.Space;
import etomica.species.Species;
import etomica.util.Arrays;

public abstract class PotentialMasterNbr extends PotentialMaster implements AtomTypeAgentManager.AgentSource {

    protected PotentialMasterNbr(Space space, PhaseAgentSource phaseAgentSource, 
            PhaseAgentManager phaseAgentManager) {
        super(space);
        this.phaseAgentSource = phaseAgentSource;
        this.phaseAgentManager = phaseAgentManager;
        rangedAgentManager = new AtomTypeAgentManager(this, null);
        intraAgentManager = new AtomTypeAgentManager(this, null);
        initialized = false;
    }
    
    public PotentialGroup makePotentialGroup(int nBody) {
        return new PotentialGroupNbr(nBody, space);
    }
    
    /**
     * Initializes the AtomType-based potential arrays and adds us as a 
     * SimulationListener so we'll know when AtomTypes are added or removed
     * (which only happens when Species are added or removed).
     */
    protected void init() {
        rangedAgentManager.setRoot(getSimulation().getSpeciesRoot());
        intraAgentManager.setRoot(getSimulation().getSpeciesRoot());
        rangedPotentialIterator = rangedAgentManager.makeIterator();
        intraPotentialIterator = intraAgentManager.makeIterator();
        initialized = true;
    }
    
    public void addPotential(Potential potential, Species[] species) {
        if (!initialized) {
            init();
        }
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
             addRangedPotentialForTypes(potential, atomTypes);
        }
    }

    public void potentialAddedNotify(Potential subPotential, PotentialGroup pGroup) {
        if (!initialized) {
            init();
        }
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
                ((PotentialArray)intraAgentManager.getAgent(parentType[0])).addPotential(pGroup);
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
        addRangedPotentialForTypes(subPotential, atomTypes);
    }

    protected abstract void addRangedPotentialForTypes(Potential subPotential, AtomType[] atomTypes);
    
    protected void addRangedPotential(Potential potential, AtomType atomType) {
        if (!initialized) {
            init();
        }
        
        PotentialArray potentialAtomType = (PotentialArray)rangedAgentManager.getAgent(atomType);
        potentialAtomType.addPotential(potential);
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
            rangedPotentialIterator.reset();
            while (rangedPotentialIterator.hasNext()) {
                ((PotentialArray)rangedPotentialIterator.next()).removePotential(potential);
            }
        }
        else if (potential instanceof PotentialGroup) {
            intraPotentialIterator.reset();
            while (intraPotentialIterator.hasNext()) {
                ((PotentialArray)intraPotentialIterator.next()).removePotential(potential);
            }
        }
        allPotentials = (Potential[])Arrays.removeObject(allPotentials,potential);
    }
    
    public PotentialArray getRangedPotentials(AtomType atomType) {
        return (PotentialArray)rangedAgentManager.getAgent(atomType);
    }

    public PotentialArray getIntraPotentials(AtomType atomType) {
        return (PotentialArray)intraAgentManager.getAgent(atomType);
    }
    
    public final PhaseAgentManager getCellAgentManager() {
        // phaseAgentManager returns performing no action if the given
        // speciesRoot is the same as was set previously.  We can't set
        // this in setSimulation because the Simulation might be not
        // have its speciesRoot yet
        phaseAgentManager.setRoot(getSimulation().getSpeciesRoot());
        return phaseAgentManager;
    }
    
    public Class getAgentClass() {
        return PotentialArray.class;
    }
    
    public Object makeAgent(AtomType type) {
        return new PotentialArray();
    }
    
    public void releaseAgent(Object agent, AtomType type) {
    }

    protected AtomTypeAgentManager.AgentIterator rangedPotentialIterator;
    protected AtomTypeAgentManager.AgentIterator intraPotentialIterator;
    protected final AtomTypeAgentManager rangedAgentManager;
    protected final AtomTypeAgentManager intraAgentManager;
    protected Potential[] allPotentials = new Potential[0];
    protected PhaseAgentSource phaseAgentSource;
    private PhaseAgentManager phaseAgentManager;
    protected boolean initialized;
}
