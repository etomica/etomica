package etomica.atom;

import etomica.atom.iterator.AtomIteratorArrayListSimple;
import etomica.atom.iterator.AtomIteratorListSimple;
import etomica.phase.Phase;
import etomica.simulation.SimulationEvent;
import etomica.simulation.SimulationEventManager;
import etomica.simulation.SimulationListener;
import etomica.species.Species;
import etomica.util.Arrays;

/**
 * Parent of all SpeciesMaster instances in a simulation, and thus the root
 * of the hierarchy of atoms.  A single instance
 * of this class is constructed and held the the Simulation. This class ensures
 * that an agent of every species is added to every phase.  
 */
 
public final class SpeciesRoot extends Atom {
    
    private final AtomTypeGroup childType;//accessed by SpeciesMaster
    private Species[] speciesList = new Species[0];
    //manager and events for addition/removal of Phases
    private final SimulationEventManager eventManager; 

    public SpeciesRoot(int[] bitLength) {
        super(new AtomTypeRoot(AtomIndexManager.makeRootIndexManager(bitLength)), new NodeFactory());
        childType = new AtomTypeGroup((AtomTypeGroup)type,null);
        node.setOrdinal(1);
        eventManager = new SimulationEventManager();
    }
    
    public String signature() {return "root";}
    
    /**
     * Adds species to the list of all species in the simulation, and
     * adds new species agent to every phase currently in simulation.
     * This is called by the Species constructor.
     * 
     * @return the index assigned to the new species
     */
    //would prefer this weren't public, but must be accessed by Species
    public int addSpecies(Species species) {
        speciesList = (Species[])Arrays.addObject(speciesList,species);
        AtomArrayList speciesMasters = ((AtomTreeNodeGroup)node).childList;
        AtomIteratorArrayListSimple iterator = new AtomIteratorArrayListSimple(speciesMasters);
        iterator.reset();
        while (iterator.hasNext()) {
            species.makeAgent((SpeciesMaster)iterator.nextAtom());
        }
        return speciesList.length;
    }

    public boolean removeSpecies(Species species) {
        boolean success = false;
        for (int i=0; i<speciesList.length; i++) {
            if (speciesList[i] == species) {
                success = true;
                break;
            }
        }
        if (!success) {
            return false;
        }
        speciesList = (Species[])Arrays.removeObject(speciesList,species);
        AtomArrayList speciesMasters = ((AtomTreeNodeGroup)node).childList;
        AtomIteratorArrayListSimple iterator = new AtomIteratorArrayListSimple(speciesMasters);
        iterator.reset();
        while (iterator.hasNext()) {
            ((SpeciesMaster)iterator.nextAtom()).removeSpecies(species);
        }
        
        AtomType[] agentTypes = childType.childTypes;
        AtomType removedType = null;
        for (int i=0; i<agentTypes.length; i++) {
            if (agentTypes[i].getSpecies() == species) {
                removedType = agentTypes[i];
            }
        }
        childType.childTypes = (AtomType[])Arrays.removeObject(childType.childTypes,removedType);
        return true;
    }
    
    public Species[] getSpecies() {
        return (Species[])speciesList.clone();
    }

    //event management
    public synchronized void addListener(SimulationListener listener) {
        eventManager.addListener(listener);
    }

    public synchronized void removeListener(SimulationListener listener) {
        eventManager.removeListener(listener);
    }


    /**
     * @return Returns the childType.
     */
    public AtomTypeGroup getChildType() {
        return childType;
    }
    
    private static final class RootAtomTreeNode extends AtomTreeNodeGroup {
        
        RootAtomTreeNode(Atom atom) {
            super(atom);
            additionEvent = new SimulationEvent(this);
            additionEvent.setType(SimulationEvent.PHASE_ADDED);
            removalEvent = new SimulationEvent(this);
            removalEvent.setType(SimulationEvent.PHASE_REMOVED);
        }
        public Phase parentPhase() {
            throw new RuntimeException("Error:  Unexpected call to parentPhase in SpeciesRoot");
        }
        
        public Species parentSpecies() {
            throw new RuntimeException("Error:  Unexpected call to parentSpecies in SpeciesRoot");
        }
        public SpeciesAgent parentSpeciesAgent() {
            throw new RuntimeException("Error:  Unexpected call to parentSpeciesAgent in SpeciesRoot");
        }
        /**
        * Returns null, because a species master is not contained within a molecule.
        */
        public final Atom parentMolecule() {
            throw new RuntimeException("Error:  Unexpected call to parentMolecule in SpeciesRoot");
        }
        
        /**
         * Ends recursive chain to determine child of given node from which this
         * node is descended.  Always returns null.
         */
        public AtomTreeNode childWhereDescendedFrom(AtomTreeNode node) {
            return null;
        }
        /**
         * Returns true, because children are SpeciesMaster instances.
         */
        public final boolean childrenAreGroups() {return true;}
        
       /**
         * Checks that new atom is a SpeciesMaster instance, and adds a
         * species agent to it for every species currently in simulation.
         */
        public void addAtomNotify(Atom newAtom) {
            if(newAtom instanceof SpeciesMaster) {
                Species[] speciesList = ((SpeciesRoot)atom).getSpecies();
                for(int i=0; i<speciesList.length; i++) {
                    speciesList[i].makeAgent((SpeciesMaster)newAtom);
                }
                additionEvent.setPhase(newAtom.node.parentPhase());
                ((SpeciesRoot)atom).eventManager.fireEvent(additionEvent);
            }
 
        }

        public void removeAtomNotify(Atom oldAtom) {
            if(oldAtom instanceof SpeciesMaster) {
//                ordinalReservoir.returnOrdinal(oldAtom.node.getOrdinal());
                removalEvent.setPhase(oldAtom.node.parentPhase());
                ((SpeciesRoot)atom).eventManager.fireEvent(removalEvent);
            }
        }
        
        private final SimulationEvent additionEvent;
        private final SimulationEvent removalEvent;
    }
    
    private static final class NodeFactory implements AtomTreeNodeFactory, java.io.Serializable {
        public AtomTreeNode makeNode(Atom atom) {
            return new RootAtomTreeNode(atom);
        }
    }
}
