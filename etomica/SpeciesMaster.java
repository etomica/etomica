package etomica;

/**
 * Coordinator of all species agents in a phase.  
 *
 * @author David Kofke
 */
 
 /* History
  * 12/06/02 (DAK) added childWhereDescendedFrom method to MasterNode subclass.
  */
public final class SpeciesMaster extends Atom {
    
    private int moleculeCount;
    //manager and events for addition/removal of descendant atoms
    private final SimulationEventManager eventManager = new SimulationEventManager();
    private final PhaseEvent additionEvent = new PhaseEvent(this, PhaseEvent.ATOM_ADDED);
    private final PhaseEvent removalEvent = new PhaseEvent(this, PhaseEvent.ATOM_REMOVED);
    public int index;
    public /*final*/ AtomTreeNodeGroup node;//shadow superclass field of same name to avoid casts
    
    /**
     * List of leaf atoms in phase, suitable for iteration via an AtomIteratorList.
     */
    public final AtomList atomList = new AtomList();

    public SpeciesMaster(Phase p) {
        super(p.parentSimulation().space(), AtomType.NULL, new NodeFactory(p), 
                p.parentSimulation().getIteratorFactory().simpleSequencerFactory(), null);//parent is null
        index = p.index;
        node = (AtomTreeNodeGroup)super.node;
  //      ((MasterAtomTreeNode)node).speciesMaster = this;
    }
        
    public void addSpecies(Species species) {
        SpeciesAgent agent = species.makeAgent(this);
        agent.setNMolecules(species.getNMolecules());
        Configuration config = node.parentPhase().getConfiguration();
        config.initializeCoordinates(node.childAtomArray());

//        parentPhase.getConfiguration().initializePositions(parentPhase.makeMoleculeIterator());
    }
    
    public SpeciesAgent firstSpecies() {return (SpeciesAgent)node.firstChildAtom();}
    public SpeciesAgent lastSpecies() {return (SpeciesAgent)node.lastChildAtom();}
        
    
//    public int atomCount() {return atomList.size();}//or could use node.leafAtomCount()
    public int moleculeCount() {return moleculeCount;}
    
    public String signature() {return node.parentPhase().getName();}
    
    
    //event management
    public synchronized void addListener(PhaseListener listener) {
        eventManager.addListener(listener);
    }

    public synchronized void removeListener(PhaseListener listener) {
        eventManager.removeListener(listener);
    }
    
    //why not make inner class?
    private static final class MasterAtomTreeNode extends AtomTreeNodeGroup {
        
        MasterAtomTreeNode(Phase parentPhase, Atom atom) {
            super(atom, null);
            speciesMaster = (SpeciesMaster)atom;
            this.parentPhase = parentPhase;
            leafIterator.setLeafIterator();
            depth = 0;
            setIndex(parentPhase.index);
        }
        public Phase parentPhase() {return parentPhase;}
        
        public Species parentSpecies() {
            throw new RuntimeException("Error:  Unexpected call to parentSpecies in SpeciesMaster");
        }
        public SpeciesAgent parentSpeciesAgent() {
            throw new RuntimeException("Error:  Unexpected call to parentSpeciesAgent in SpeciesMaster");
        }
        public Simulation parentSimulation() {return parentPhase.parentSimulation();}
        /**
        * Returns null, because a species master is not contained within a molecule.
        */
        public final Atom parentMolecule() {
            throw new RuntimeException("Error:  Unexpected call to parentMolecule in SpeciesMaster");
        }
        
        /**
         * Ends recursive chain to determine child of given node from which this
         * node is descended.  Always returns null.
         */
        public AtomTreeNode childWhereDescendedFrom(AtomTreeNode node) {
            return null;
        }
        /**
         * Returns true, because children are SpeciesAgent instances.
         */
        public final boolean childrenAreGroups() {return true;}
        
        public void addAtomNotify(Atom atom) {
            if(atom.node.parentGroup() instanceof SpeciesAgent) {speciesMaster.moleculeCount++;}
            else if(atom instanceof SpeciesAgent) {speciesMaster.moleculeCount += ((SpeciesAgent)atom).moleculeCount();}
            leafAtomCount += atom.node.leafAtomCount();
            leafIterator.setBasis(atom);
            leafIterator.reset();
            while(leafIterator.hasNext()) {
                speciesMaster.atomList.add(((AtomTreeNodeLeaf)leafIterator.next().node).leafLinker);
            }
            speciesMaster.eventManager.fireEvent(speciesMaster.additionEvent.setAtom(atom));
        }

        //updating of leaf atomList may not be efficient enough for repeated use, but is probably ok
        public void removeAtomNotify(Atom atom) {
            if(atom.node.parentGroup() instanceof SpeciesAgent) {speciesMaster.moleculeCount--;}
            else if(atom instanceof SpeciesAgent) {speciesMaster.moleculeCount -= ((SpeciesAgent)atom).moleculeCount();}
            leafAtomCount -= atom.node.leafAtomCount();
            leafIterator.setBasis(atom);
            leafIterator.reset();
            while(leafIterator.hasNext()) {
                speciesMaster.atomList.remove(((AtomTreeNodeLeaf)leafIterator.next().node).leafLinker);
            }
            speciesMaster.eventManager.fireEvent(speciesMaster.removalEvent.setAtom(atom));
        }
        
        private final Phase parentPhase;
        private final SpeciesMaster speciesMaster;
        private final AtomIteratorTree leafIterator = new AtomIteratorTree();
    } //end of MasterAtomTreeNode           
    
    private static final class NodeFactory implements AtomTreeNode.Factory {
        Phase phase;
        NodeFactory(Phase p) {
            phase = p;
        }
        public AtomTreeNode makeNode(Atom atom, AtomTreeNodeGroup dummy) {
            return new MasterAtomTreeNode(phase, atom);
        }
    }


    /**
     * non-graphic main method to test handling of leaf atom list.
     */
    public static void main(String args[]) {
        
        Simulation sim = new Simulation();
        Simulation.instance = sim;
        Species species2 = new SpeciesSpheresMono();
        Species species1 = new SpeciesSpheres(3,3);
        Species species0 = new SpeciesSpheres(3,2);
        species0.setNMolecules(4);
        species1.setNMolecules(2);
        species2.setNMolecules(2);
        Phase phase = new Phase();
        sim.elementCoordinator.go();
        
        AtomIteratorList listIterator = new AtomIteratorList();
        listIterator.setBasis(phase.speciesMaster.atomList);
        listIterator.reset();
        while(listIterator.hasNext()) System.out.println(listIterator.next().toString());
        System.out.println();
    }//end of main

}//end of SpeciesMaster