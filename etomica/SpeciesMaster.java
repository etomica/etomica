package etomica;

/**
 * Coordinator of all species agents in a phase.  
 *
 * @author David Kofke
 */
public final class SpeciesMaster extends AtomGroup {
    
    private final Phase parentPhase;
    private int moleculeCount;
    //manager and events for addition/removal of descendant atoms
    private final SimulationEventManager eventManager = new SimulationEventManager();
    private final PhaseEvent additionEvent = new PhaseEvent(this, PhaseEvent.ATOM_ADDED);
    private final PhaseEvent removalEvent = new PhaseEvent(this, PhaseEvent.ATOM_REMOVED);
    public final int index;
    
    public SpeciesMaster(Phase p) {
        super(p.parentSimulation().space(), AtomType.NULL);
        parentPhase = p;
        index = p.index;
    }
        
    public void addSpecies(Species species) {
        SpeciesAgent agent = species.makeAgent(this);
        addAtom(agent);
        parentPhase.getConfiguration().initializeCoordinates(childAtomArray());

//        parentPhase.getConfiguration().initializePositions(parentPhase.makeMoleculeIterator());
    }
    
    public int depth() {return 0;}
    
    public SpeciesAgent firstSpecies() {return (SpeciesAgent)firstChildAtom();}
    public SpeciesAgent lastSpecies() {return (SpeciesAgent)lastChildAtom();}
        
    public Phase parentPhase() {return parentPhase;}
    public Species parentSpecies() {return null;} //never called
    public SpeciesAgent parentSpeciesAgent() {return null;}//never called
    public Simulation parentSimulation() {return parentPhase.parentSimulation();}
    
    protected void addAtomNotify(Atom atom) {
        if(atom.parentGroup() instanceof SpeciesAgent) {moleculeCount++;}
        else if(atom instanceof SpeciesAgent) {moleculeCount += ((SpeciesAgent)atom).moleculeCount();}
        leafAtomCount += atom.leafAtomCount();
        eventManager.fireEvent(additionEvent.setAtom(atom));
    }

    protected void removeAtomNotify(Atom atom) {
        if(atom.parentGroup() instanceof SpeciesAgent) {moleculeCount--;}
        else if(atom instanceof SpeciesAgent) {moleculeCount -= ((SpeciesAgent)atom).moleculeCount();}
        leafAtomCount -= atom.leafAtomCount();
        eventManager.fireEvent(removalEvent.setAtom(atom));
    }
    
    public int atomCount() {return leafAtomCount;}
    public int moleculeCount() {return moleculeCount;}
    
    public String signature() {return parentPhase.getName();}
    
    
    //event management
    public synchronized void addListener(PhaseEventListener listener) {
        eventManager.addListener(listener);
    }

    public synchronized void removeListener(PhaseEventListener listener) {
        eventManager.removeListener(listener);
    }

}//end of SpeciesMaster