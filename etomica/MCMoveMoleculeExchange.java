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
    private PotentialMaster.Agent phasePotential;
    private final IteratorDirective iteratorDirective = new IteratorDirective(IteratorDirective.BOTH);
    private final PotentialCalculation.EnergySum energy = new PotentialCalculation.EnergySum();

    public MCMoveMoleculeExchange(IntegratorMC parent) {
        super();
        parentIntegrator = parent;
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
        super.setPhase(p);
        firstPhase = p[0];
        secondPhase = p[1];
    }
        
    public void thisTrial() {
        Phase iPhase, dPhase;
        if(Simulation.random.nextDouble() < 0.5) {
            iPhase = firstPhase;
            dPhase = secondPhase;
        }
        else {
            iPhase = secondPhase;
            dPhase = firstPhase;
        }
        if(dPhase.moleculeCount() == 0) {return;} //no molecules to delete; trial is over
        
        Atom a = dPhase.randomMolecule();  //select random molecule to delete
        Species species = a.parentSpecies();
        
        SpeciesAgent iSpecies = species.getAgent(iPhase);  //insertion-phase speciesAgent
        SpeciesAgent dSpecies = species.getAgent(dPhase);  //deletion-phase species Agent
        
        double uOld = dPhase.potential.calculate(iteratorDirective.set(a), energy.reset()).sum();

        a.coord.displaceTo(iPhase.randomPosition());         //place at random in insertion phase
        dPhase.removeMolecule(a,dSpecies);
        iPhase.addMolecule(a,iSpecies);//this must be done after the displaceTo call, because addMolecule may impose PBC, which would cause to forget the original position
        iPhase.iteratorFactory().moveNotify(a);
        double uNew = iPhase.potential.calculate(iteratorDirective, energy.reset()).sum();
        if(uNew == Double.MAX_VALUE) {  //overlap; reject
            a.coord.replace();                //put it back 
            iPhase.removeMolecule(a,iSpecies);
            dPhase.addMolecule(a,dSpecies);
            return;        
        }
        
        double bFactor = (dSpecies.moleculeCount()+1)/dPhase.volume()  //acceptance probability
                         * iPhase.volume()/(iSpecies.moleculeCount())                   //note that dSpecies.nMolecules has been decremented
                         * Math.exp(-(uNew-uOld)/parentIntegrator.temperature);    //and iSpecies.nMolecules has been incremented
        if(bFactor > 1.0 || bFactor > Simulation.random.nextDouble()) {  //accept
            nAccept++;
        }
        else {              //reject
            a.coord.replace();    //put it back
            iPhase.removeMolecule(a,iSpecies);
            dPhase.addMolecule(a,dSpecies); //place molecule in phase
        }
    }//end of thisTrial
}//end of MCMoveMoleculeExchange