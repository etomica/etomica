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
    
    protected Integrator integrator;
    
    public SpeciesAgent(Species s, int nMolecules) {
        super(s.parentSimulation().space(), AtomType.NULL, new NodeFactory(s));
        factory = s.moleculeFactory();
        for(int i=0; i<nMolecules; i++) {
            node.addAtom(factory.makeAtom());
        }
        node.setDepth(1);
    }
        
    public final AtomFactory moleculeFactory() {return factory;}
    
    
    public SpeciesAgent nextSpecies() {return (SpeciesAgent)nextAtom();}
    public int moleculeCount() {return node.childAtomCount();}
    public Atom firstMolecule() {return node.firstChildAtom();}
    public Atom lastMolecule() {return node.lastChildAtom();}
    public Atom randomMolecule() {return node.randomAtom();}
    
//    public final int depth() {return 1;}
        
    public Atom addNewAtom() {
        if(!resizable) return null; //should define an exeception 
        Atom aNew = moleculeFactory().makeAtom();
        node.addAtom(aNew);
        return aNew;
    }
    
    /**
     * Performs the given action on all children (molecules) of this species agent.
     */
    public void allMolecules(AtomAction action) {
        Atom last = node.lastChildAtom();
        for(Atom a=node.firstChildAtom(); a!=null; a=a.nextAtom()) {
            action.actionPerformed(a);
            if(a == last) break;
        }
    }
    
    /**
     * Performs the given action on all (leaf) atoms of this species agent.
     */
    public void allAtoms(AtomAction action) {
        Atom last = node.lastLeafAtom();
        for(Atom a=node.firstLeafAtom(); a!=null; a=a.nextAtom()) {
            action.actionPerformed(a);
            if(a == last) break;
        }
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
        if(!resizable) return;
        boolean wasPaused = pauseIntegrator();
        if(forceRebuild) node.removeAll();
        
        if(n <= 0) node.removeAll();
        else if(n > node.childAtomCount()) {
            for(int i=node.childAtomCount(); i<n; i++) node.addAtom(factory.makeAtom());
        }
        else if(n < node.childAtomCount()) {
            for(int i=node.childAtomCount(); i>n; i--) node.removeAtom(node.lastChildAtom());
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
        AgentAtomTreeNode(Species parentSpecies, Atom atom) {
            super(atom);
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
        public AtomTreeNode makeNode(Atom atom) {
            return new AgentAtomTreeNode(species, atom);
        }
    }

} //end of SpeciesAgent
