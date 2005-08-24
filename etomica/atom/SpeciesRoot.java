package etomica.atom;

import etomica.atom.iterator.AtomIteratorListSimple;
import etomica.phase.Phase;
import etomica.phase.PhaseEvent;
import etomica.phase.PhaseListener;
import etomica.simulation.SimulationEventManager;
import etomica.species.Species;
import etomica.util.Arrays;

/**
 * Parent of all SpeciesMaster instances in a simulation, and thus the root
 * of the hierarcy of atoms.  A single instance
 * of this class is constructed and held the the Simulation. This class ensures
 * that an agent of every species is added to every phase.  
 */
 
public final class SpeciesRoot extends Atom {
    
    private final AtomTypeGroup childType;//accessed by SpeciesMaster
    private Species[] speciesList = new Species[0];
    //manager and events for addition/removal of descendant atoms
    private final SimulationEventManager eventManager = new SimulationEventManager();
    private final PhaseEvent additionEvent = new PhaseEvent(this,
            PhaseEvent.ATOM_ADDED);
    private final PhaseEvent removalEvent = new PhaseEvent(this,
            PhaseEvent.ATOM_REMOVED);

    public SpeciesRoot(int[] bitLength) {
        super(null, new AtomTypeRoot(AtomIndexManager.makeRootIndexManager(bitLength)), new NodeFactory(), AtomLinker.FACTORY);
        childType = new AtomTypeGroup((AtomTypeGroup)type,null);
        node.setOrdinal(1);
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
        AtomList speciesMasters = ((AtomTreeNodeGroup)node).childList;
        AtomIteratorListSimple iterator = new AtomIteratorListSimple(speciesMasters);
        iterator.reset();
        while (iterator.hasNext()) {
            species.makeAgent((SpeciesMaster)iterator.nextAtom());
        }
        return speciesList.length;
    }

    public void removeSpecies(Species species) {
        speciesList = (Species[])Arrays.removeObject(speciesList,species);
        AtomList speciesMasters = ((AtomTreeNodeGroup)node).childList;
        AtomIteratorListSimple iterator = new AtomIteratorListSimple(speciesMasters);
        iterator.reset();
        while (iterator.hasNext()) {
//            ((SpeciesMaster)iterator.nextAtom()).removeSpecies(species);
        }
    }

    //event management
    public synchronized void addListener(PhaseListener listener) {
        eventManager.addListener(listener);
    }

    public synchronized void removeListener(PhaseListener listener) {
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
                Species[] speciesList = ((SpeciesRoot)atom).speciesList;
                for(int i=0; i<speciesList.length; i++) {
                    speciesList[i].makeAgent((SpeciesMaster)newAtom);
                }
            }
            ((SpeciesRoot)atom).eventManager.fireEvent(((SpeciesRoot)atom).additionEvent
                    .setAtom(atom));
 
        }

        public void removeAtomNotify(Atom atom) {
            ((SpeciesRoot)atom).eventManager.fireEvent(((SpeciesRoot)atom).removalEvent
                    .setAtom(atom));
        }
    }
    
    private static final class NodeFactory implements AtomTreeNodeFactory, java.io.Serializable {
        public AtomTreeNode makeNode(Atom atom) {
            return new RootAtomTreeNode(atom);
        }
    }
}
