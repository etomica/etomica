package etomica;

/**
 * The SpeciesAgent is a representative of the species in each phase.
 * The agent handles addition, deletion, link-list ordering, counting, etc. of 
 * molecules in a phase.  Each phase has an agent from every species instance.
 * 
 * @author David Kofke
 */

public class SpeciesAgent extends AtomGroup {

    private final Species parentSpecies;
    protected final AtomFactory factory;
    
    protected Integrator integrator;
    
    public SpeciesAgent(Species s, SpeciesMaster parent, int nMolecules) {
        super(parent, AtomType.NULL);
        parentSpecies = s;
        factory = s.moleculeFactory();
        for(int i=0; i<nMolecules; i++) {
            addAtom(factory.makeAtom());
        }
    }
        
    public AtomFactory atomFactory() {return factory;}
    
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
    public final SpeciesAgent parentSpeciesAgent() {return this;}
    
    public SpeciesAgent nextSpecies() {return (SpeciesAgent)nextAtom();}
    public int moleculeCount() {return childCount();}
    public Atom firstMolecule() {return firstChild();}
    public Atom lastMolecule() {return lastChild();}
    
    public final int depth() {return 1;}
        
    public Atom addNewAtom() {
        if(!resizable) return null; //should define an exeception 
        Atom aNew = atomFactory().makeAtom();
        addAtom(aNew);
        return aNew;
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
        if(forceRebuild) removeAll();
        
        if(n <= 0) removeAll();
        else if(n > childCount) {
            for(int i=childCount; i<n; i++) addAtom(factory.makeAtom());
        }
        else if(n < childCount) {
            for(int i=childCount; i>n; i--) removeAtom(lastChild());
        }
        
        //reconsider this
        parentPhase().configuration.initializeCoordinates(this);
        parentPhase().iteratorFactory().reset();
        
        unpauseIntegrator(wasPaused);
    }
    
    private boolean pauseIntegrator() {
        Phase phase = parentPhase();
        integrator = (phase != null) ? phase.integrator() : null;
        boolean wasPaused = true;
        if(integrator != null) {
            wasPaused = integrator.isPaused();//record pause state of integrator
            if(!wasPaused) {
                integrator.pause();
                while(!integrator.isPaused()) {}
            }
        }
        return wasPaused;
    }
    
    private void unpauseIntegrator(boolean wasPaused) {
        if(integrator != null) {
            if(integrator.isInitialized()) integrator.initialize();//reinitialize only if initialized already
            if(!wasPaused) integrator.unPause();//resume if was not paused originally
        }
    }
              
} //end of SpeciesAgent
