package etomica;

/**
 * Performs a trial that results in the exchange of a molecule from one phase to another.
 * Primary use is as an elementary move in a Gibbs ensemble simulation
 *
 * @author David Kofke
 */
public final class MCMoveMoleculeExchange extends MCMove {
    
    private Phase firstPhase;
    private Phase secondPhase;
    private final double ROOT;
    private final IteratorDirective iteratorDirective = new IteratorDirective(IteratorDirective.BOTH);
    private final AtomIteratorSinglet affectedAtomIterator = new AtomIteratorSinglet();

    private transient Atom molecule;
    private transient Phase iPhase, dPhase;
    private transient SpeciesAgent iSpecies, dSpecies;
    private transient double uOld;

    public MCMoveMoleculeExchange(IntegratorMC parent) {
        super(parent);
        ROOT = 1.0/(double)parentIntegrator.parentSimulation().space().D();
        setTunable(false);
        setPerParticleFrequency(true);
    }
    
    /**
     * Overrides superclass method so that it performs no action.
     * Must set using method that takes an array of phases.
     */
    public void setPhase(Phase p) {}

    public void setPhase(Phase[] p) {
        if(p == null || p.length == 0) return;
        super.setPhase(p);
        firstPhase = p[0];
        if(p.length < 2) return;
        secondPhase = p[1];
    }
        
    public boolean doTrial() {
        if(Simulation.random.nextDouble() < 0.5) {
            iPhase = firstPhase;
            dPhase = secondPhase;
        }
        else {
            iPhase = secondPhase;
            dPhase = firstPhase;
        }
        if(dPhase.moleculeCount() == 0) {return false;} //no molecules to delete; trial is over
        
        molecule = dPhase.randomMolecule();  //select random molecule to delete
        Species species = molecule.node.parentSpecies();
        
        iSpecies = species.getAgent(iPhase);  //insertion-phase speciesAgent
        dSpecies = species.getAgent(dPhase);  //deletion-phase species Agent
        
        uOld = potential.set(dPhase).calculate(iteratorDirective.set(molecule), energy.reset()).sum();

        molecule.coord.displaceTo(iPhase.randomPosition());         //place at random in insertion phase
        molecule.node.setParent(iSpecies);
        return true;
    }//end of doTrial
    
    public double lnTrialRatio() {
        //note that dSpecies.nMolecules has been decremented
        //and iSpecies.nMolecules has been incremented
        return Math.log( (dSpecies.moleculeCount()+1)/dPhase.volume()
                         * iPhase.volume()/iSpecies.moleculeCount() ); 
    }
    
    public double lnProbabilityRatio() {
        double uNew = potential.set(iPhase).calculate(iteratorDirective, energy.reset()).sum();
        return -(uNew - uOld)/parentIntegrator.temperature;
    }
    
    public void acceptNotify() {  /* do nothing */}
    
    public void rejectNotify() {
        molecule.coord.replace();
        molecule.node.setParent(dSpecies);
    }

    public final AtomIterator affectedAtoms(Phase phase) {
        if(this.firstPhase != phase && this.secondPhase != phase) return AtomIterator.NULL;
        affectedAtomIterator.setBasis(molecule);
        affectedAtomIterator.reset();
        return affectedAtomIterator;
    }

}//end of MCMoveMoleculeExchange