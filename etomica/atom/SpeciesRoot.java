package etomica.atom;

import etomica.atom.iterator.AtomIteratorArrayListSimple;
import etomica.phase.Phase;
import etomica.simulation.SimulationEventManager;
import etomica.simulation.SimulationPhaseAddedEvent;
import etomica.simulation.SimulationPhaseRemovedEvent;
import etomica.simulation.SimulationSpeciesAddedEvent;
import etomica.simulation.SimulationSpeciesRemovedEvent;
import etomica.species.Species;
import etomica.util.Arrays;

/**
 * Parent of all SpeciesMaster instances in a simulation, and thus the root
 * of the hierarchy of atoms.  A single instance
 * of this class is constructed and held the the Simulation. This class ensures
 * that an agent of every species is added to every phase.  
 */
 
public final class SpeciesRoot extends Atom {
    
    public SpeciesRoot(int[] bitLength) {
        super(new AtomTypeRoot(AtomAddressManager.makeRootIndexManager(bitLength)), new NodeFactory());
        speciesMasterType = new AtomTypeGroup((AtomPositionDefinition)null);
        speciesMasterType.setParentType((AtomTypeGroup)type);
        node.setIndex(0);
        eventManager = new SimulationEventManager();
        ((AtomTypeRoot)type).setEventManager(eventManager);
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
        for (int i=0; i<speciesList.length; i++) {
            if (speciesList[i] == species) {
                throw new IllegalArgumentException("Species already exists");
            }
        }
        speciesList = (Species[])Arrays.addObject(speciesList,species);
        AtomArrayList speciesMasters = ((AtomTreeNodeGroup)node).childList;
        AtomIteratorArrayListSimple iterator = new AtomIteratorArrayListSimple(speciesMasters);
        iterator.reset();
        while (iterator.hasNext()) {
            species.makeAgent((SpeciesMaster)iterator.nextAtom());
        }

        eventManager.fireEvent(new SimulationSpeciesAddedEvent(species));
        
        return speciesList.length-1;
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

        eventManager.fireEvent(new SimulationSpeciesRemovedEvent(species));
        
        speciesList = (Species[])Arrays.removeObject(speciesList,species);
        AtomArrayList speciesMasters = ((AtomTreeNodeGroup)node).childList;
        AtomIteratorArrayListSimple iterator = new AtomIteratorArrayListSimple(speciesMasters);
        iterator.reset();
        while (iterator.hasNext()) {
            ((SpeciesMaster)iterator.nextAtom()).removeSpecies(species);
        }
        
        ((AtomTypeRoot)type).removeSpecies(species);
        return true;
    }
    
    public Species[] getSpecies() {
        return (Species[])speciesList.clone();
    }

    public SimulationEventManager getEventManager() {
        return eventManager;
    }
    
    /**
     * @return Returns the AtomType of the SpeciesMaster
     */
    public AtomTypeGroup getSpeciesMasterType() {
        return speciesMasterType;
    }
    
    private static final long serialVersionUID = 1L;
    private final AtomTypeGroup speciesMasterType;//accessed by SpeciesMaster
    private Species[] speciesList = new Species[0];
    //manager and events for addition/removal of Phases
    protected final SimulationEventManager eventManager;

    private static final class RootAtomTreeNode extends AtomTreeNodeGroup {
        
        RootAtomTreeNode(Atom atom) {
            super(atom);
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
                ((SpeciesRoot)atom).eventManager.fireEvent(new SimulationPhaseAddedEvent(newAtom.getNode().parentPhase()));
            }
 
        }

        public void removeAtomNotify(Atom oldAtom) {
            if(oldAtom instanceof SpeciesMaster) {
//                ordinalReservoir.returnOrdinal(oldAtom.node.getOrdinal());
                ((SpeciesRoot)atom).eventManager.fireEvent(new SimulationPhaseRemovedEvent(oldAtom.getNode().parentPhase()));
            }
        }
        
        private static final long serialVersionUID = 1L;
    }
    
    protected static final class NodeFactory implements AtomTreeNodeFactory, java.io.Serializable {
        public AtomTreeNode makeNode(Atom atom) {
            return new RootAtomTreeNode(atom);
        }
        private static final long serialVersionUID = 1L;
    }
}
