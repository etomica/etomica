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
    private final AtomIteratorSequential affectedAtomIterator = new AtomIteratorSequential(true);
    private Atom molecule;

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
        
    public boolean thisTrial() {
        Phase iPhase, dPhase;
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
        
        SpeciesAgent iSpecies = species.getAgent(iPhase);  //insertion-phase speciesAgent
        SpeciesAgent dSpecies = species.getAgent(dPhase);  //deletion-phase species Agent
        
        double uOld = potential.set(dPhase).calculate(iteratorDirective.set(molecule), energy.reset()).sum();

        molecule.coord.displaceTo(iPhase.randomPosition());         //place at random in insertion phase
        molecule.node.setParent(iSpecies);
 //       dPhase.removeMolecule(molecule,dSpecies);
 //       iPhase.addMolecule(molecule,iSpecies);//this must be done after the displaceTo call, because addMolecule may impose PBC, which would cause to forget the original position
        double uNew = potential.set(iPhase).calculate(iteratorDirective, energy.reset()).sum();
        if(uNew == Double.MAX_VALUE) {  //overlap; reject
            molecule.coord.replace();                //put it back 
            molecule.node.setParent(dSpecies);
        //    iPhase.removeMolecule(molecule,iSpecies);
        //    dPhase.addMolecule(molecule,dSpecies);
            return false;        
        }
        
        double bFactor = (dSpecies.moleculeCount()+1)/dPhase.volume()  //acceptance probability
                         * iPhase.volume()/(iSpecies.moleculeCount())                   //note that dSpecies.nMolecules has been decremented
                         * Math.exp(-(uNew-uOld)/parentIntegrator.temperature);    //and iSpecies.nMolecules has been incremented
        if(bFactor > 1.0 || bFactor > Simulation.random.nextDouble()) {  //accept
            return true;
        }
        else {              //reject
            molecule.coord.replace();    //put it back
            molecule.node.setParent(dSpecies);
        //    iPhase.removeMolecule(molecule,iSpecies);
        //    dPhase.addMolecule(molecule,dSpecies); //place molecule in phase
            return false;
        }
    }//end of thisTrial
    
    public final AtomIterator affectedAtoms() {
        affectedAtomIterator.setBasis(molecule);
        return affectedAtomIterator;
    }

}//end of MCMoveMoleculeExchange