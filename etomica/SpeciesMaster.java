package etomica;

/**
 * Coordinator of all species agents in a phase.  
 *
 * @author David Kofke
 */
public final class SpeciesMaster extends AtomGroup {
    
    private int moleculeCount;
    //manager and events for addition/removal of descendant atoms
    private final SimulationEventManager eventManager = new SimulationEventManager();
    private final PhaseEvent additionEvent = new PhaseEvent(this, PhaseEvent.ATOM_ADDED);
    private final PhaseEvent removalEvent = new PhaseEvent(this, PhaseEvent.ATOM_REMOVED);
    public int index;
    
    public SpeciesMaster(Phase p) {
        super(p.parentSimulation().space(), AtomType.NULL, new AtomTreeNode(p));
        index = p.index;
        ((AtomTreeNode)node).speciesMaster = this;
    }
        
    public void addSpecies(Species species) {
        SpeciesAgent agent = species.makeAgent(this);
        node.addAtom(agent);
        Configuration config = node.parentPhase().getConfiguration();
        config.initializeCoordinates(node.childAtomArray());

//        parentPhase.getConfiguration().initializePositions(parentPhase.makeMoleculeIterator());
    }
    
    public SpeciesAgent firstSpecies() {return (SpeciesAgent)node.firstChildAtom();}
    public SpeciesAgent lastSpecies() {return (SpeciesAgent)node.lastChildAtom();}
        
    
    public int atomCount() {return node.leafAtomCount();}
    public int moleculeCount() {return moleculeCount;}
    
    public String signature() {return node.parentPhase().getName();}
    
    
    //event management
    public synchronized void addListener(PhaseListener listener) {
        eventManager.addListener(listener);
    }

    public synchronized void removeListener(PhaseListener listener) {
        eventManager.removeListener(listener);
    }
    
    private static final class AtomTreeNode extends AtomTreeNodeGroup {
        
        private final Phase parentPhase;
        SpeciesMaster speciesMaster;
        AtomTreeNode(Phase parentPhase) {
            this.parentPhase = parentPhase;
            depth = 0;
        }
        public Phase parentPhase() {return parentPhase;}
        public Species parentSpecies() {return null;} //never called
        public SpeciesAgent parentSpeciesAgent() {return null;}//never called
        public Simulation parentSimulation() {return parentPhase.parentSimulation();}
        /**
        * Returns null, since a species master is not contained within a molecule.
        */
        public final Atom parentMolecule() {return null;}
        
        public void addAtomNotify(Atom atom) {
            if(atom.node.parentGroup() instanceof SpeciesAgent) {speciesMaster.moleculeCount++;}
            else if(atom instanceof SpeciesAgent) {speciesMaster.moleculeCount += ((SpeciesAgent)atom).moleculeCount();}
            leafAtomCount += atom.node.leafAtomCount();
            speciesMaster.eventManager.fireEvent(speciesMaster.additionEvent.setAtom(atom));
        }

        public void removeAtomNotify(Atom atom) {
            if(atom.node.parentGroup() instanceof SpeciesAgent) {speciesMaster.moleculeCount--;}
            else if(atom instanceof SpeciesAgent) {speciesMaster.moleculeCount -= ((SpeciesAgent)atom).moleculeCount();}
            leafAtomCount -= atom.node.leafAtomCount();
            speciesMaster.eventManager.fireEvent(speciesMaster.removalEvent.setAtom(atom));
        }
    }            

}//end of SpeciesMaster