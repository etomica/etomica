package etomica;
import etomica.units.Dimension;

/**
 * The SpeciesAgent is a representative of the species in each phase.
 * The agent handles addition, deletion, link-list ordering, counting, etc. of 
 * molecules in a phase.  Each phase has an agent from every species instance.
 * 
 * @author David Kofke
 */

public final class SpeciesAgent extends AtomGroup {

    protected final AtomFactory factory;
//    public final AtomTreeNodeGroup node;//shadow superclass field of same name to avoid casts
    private final AtomIteratorTree leafIterator = new AtomIteratorTree(this);
    private final AtomIteratorList moleculeIterator = new AtomIteratorList(((AtomTreeNodeGroup)this.node).childList);
    
    protected Integrator integrator;
    
    public SpeciesAgent(Species s, int nMolecules, AtomTreeNodeGroup parent) {
        super(s.parentSimulation().space(), AtomType.NULL, new NodeFactory(s), parent);
        factory = s.moleculeFactory();
//        node = (AtomTreeNodeGroup)super.node;
    }
        
    public final AtomFactory moleculeFactory() {return factory;}
    
    
    public SpeciesAgent nextSpecies() {return (SpeciesAgent)seq.next.atom;}
    public int moleculeCount() {return ((AtomTreeNodeGroup)node).childAtomCount();}
    public Atom firstMolecule() {return ((AtomTreeNodeGroup)node).firstChildAtom();}
    public Atom lastMolecule() {return ((AtomTreeNodeGroup)node).lastChildAtom();}
    public Atom randomMolecule() {return ((AtomTreeNodeGroup)node).randomAtom();}
            
    public Atom addNewAtom() {
        if(!resizable) return null; //should define an exeception 
        Atom aNew = moleculeFactory().makeAtom((AtomTreeNodeGroup)this.node);
        return aNew;
    }
    
    /**
     * Performs the given action on all children (molecules) of this species agent.
     * Method provides the iterator that performs the action, so it must be
     * synchronized to prevent multiple processes from using the single iterator.
     */
    public synchronized void allMolecules(AtomAction action) {
        moleculeIterator.allAtoms(action);
    }
    
    /**
     * Performs the given action on all (leaf) atoms of this species agent.
     * Method provides the iterator that performs the action, so it must be
     * synchronized to prevent multiple processes from using the single iterator.
     */
    public synchronized void allAtoms(AtomAction action) {
        leafIterator.allAtoms(action);
    }
    
    /**
    * Sets the number of molecules for this species.  Makes the given number
    * of new molecules, linked-list orders and initializes them.
    * Any previously existing molecules for this species in this phase are abandoned
    * Any links to molecules of next or previous species are maintained.
    * Takes no action at all if the new number of molecules equals the existing number
    *
    * @param n  the new number of molecules for this species
    * @see #makeMolecule
    * @see #deleteMolecule
    * @see #addMolecule
    */
    public void setNMolecules(int n) {
        setNMolecules(n, false);
    }
        
    /**
     * Same as setNMolecules, but takes a boolean argument that can indicate that all
     * new molecules should be made.
     */
    public void setNMolecules(int n, boolean forceRebuild) {
        AtomTreeNodeGroup treeNode = (AtomTreeNodeGroup)node;
        if(!resizable) return;
        boolean wasPaused = pauseIntegrator();
        if(forceRebuild) while(treeNode.childList.size() > 0) treeNode.lastChildAtom().sendToReservoir();
        
        else if(n > treeNode.childAtomCount()) {
            for(int i=treeNode.childAtomCount(); i<n; i++) addNewAtom();
        }
        else if(n < treeNode.childAtomCount()) {
            if(n < 0) n = 0;
            for(int i=treeNode.childAtomCount(); i>n; i--) treeNode.lastChildAtom().sendToReservoir();
        }
        
        //reconsider this
        node.parentPhase().configuration.initializeCoordinates(this);
        
        unpauseIntegrator(wasPaused);
    }
    
    public int getNMolecules() {return moleculeCount();}
    public Dimension getNMoleculesDimension() {return Dimension.QUANTITY;}

    private boolean pauseIntegrator() {
        Phase phase = node.parentPhase();
        integrator = (phase != null) ? phase.integrator() : null;
        boolean wasPaused = true;
        if(integrator != null) {
            wasPaused = integrator.isPaused();//record pause state of integrator
            if(!wasPaused) integrator.pause();//and waits until integrator puts pause in effect
        }
        return wasPaused;
    }
    
    private void unpauseIntegrator(boolean wasPaused) {
        if(integrator != null) {
            if(integrator.isInitialized()) integrator.initialize();//reinitialize only if initialized already
            if(!wasPaused) integrator.unPause();//resume if was not paused originally
        }
    }
    
    /**
     * Special AtomTreeNode class for SpeciesAgent.
     */
    private static final class AgentAtomTreeNode extends AtomTreeNodeGroup {
        
        private final Species parentSpecies;
        AgentAtomTreeNode(Species parentSpecies, Atom atom, AtomTreeNodeGroup speciesMasterNode) {
            super(atom, speciesMasterNode);
            this.parentSpecies = parentSpecies;
            depth = 1;
        }

        /**
        * Overrides super class method and terminates recursive call to identify
        * a constituent atom's species.
        */
        public final Species parentSpecies() {return parentSpecies;}
        /**
        * Overrides super class method and terminates recursive call to identify
        * a constituent atom's simulation.
        */
        public Simulation parentSimulation() {return parentSpecies.parentSimulation();}    
        /**
        * Overrides parent class method and terminates recursive call to identify this
        * as a constituent atom's species agent.
        */
        public final SpeciesAgent parentSpeciesAgent() {return (SpeciesAgent)this.atom;}
        
        /**
        * Returns null, since a species agent is not contained within a molecule.
        */
        public final Atom parentMolecule() {return null;}
    }
    
    private static final class NodeFactory implements AtomTreeNode.Factory {
        Species species;
        NodeFactory(Species s) {
            species = s;
        }
        public AtomTreeNode makeNode(Atom atom, AtomTreeNodeGroup speciesMasterNode) {
            return new AgentAtomTreeNode(species, atom, speciesMasterNode);
        }
    }

} //end of SpeciesAgent
