package etomica;

/**
 * Elementary Monte Carlo move in which a molecule of a specified species is
 * inserted into or removed from a phase.
 *
 * @author David Kofke
 */
public class MCMoveInsertDelete extends MCMove {
    
    //chemical potential
    private double mu;
    
    //directive must specify "BOTH" to get energy with all atom pairs
    private final IteratorDirective iteratorDirective = new IteratorDirective(IteratorDirective.BOTH);
    private Species species;
    private SpeciesAgent speciesAgent;
    private final AtomIteratorSinglet affectedAtomIterator = new AtomIteratorSinglet();
    private Atom testMolecule;
    private double hOld;
    private boolean insert;
    private AtomReservoir reservoir;

    public MCMoveInsertDelete(IntegratorMC parent) {
        super(parent);
        setStepSizeMax(1.0);
        setStepSizeMin(0.0);
        setStepSize(0.10);
        setMu(0.0);
        setTunable(false);
    }
    
//perhaps should have a way to ensure that two instances of this class aren't assigned the same species
    public void setSpecies(Species s) {
        species = s;
        if(phase != null) speciesAgent = (SpeciesAgent)species.getAgent(phase); 
        reservoir = new AtomReservoir(s.moleculeFactory());
    }
    public Species getSpecies() {return species;}
    
    public void setPhase(Phase p) {
        super.setPhase(p);
        if(species != null) speciesAgent = (SpeciesAgent)species.getAgent(phase); 
    }
    
    /**
     * Chooses and performs with equal probability an elementary molecule insertion
     * or deletion.
     */
    public boolean doTrial() {
        insert = Simulation.random.nextDouble() < 0.5;
        if(insert) {
            hOld = 0.0;
            testMolecule = species.moleculeFactory().makeAtom((AtomTreeNodeGroup)speciesAgent.node);
            testMolecule.coord.translateTo(phase.randomPosition());
        } else {//delete
            if(speciesAgent.moleculeCount() == 0) return false;
            testMolecule = speciesAgent.randomMolecule();
            hOld = -mu + potential.set(phase).calculate(iteratorDirective.set(testMolecule), energy.reset()).sum();
            reservoir.addAtom(testMolecule);
        }    
        return true;
    }//end of doTrial
    
    public double lnTrialRatio() {//note that moleculeCount() gives the number of molecules after the trial is attempted
        return insert ? Math.log(phase.volume()/speciesAgent.moleculeCount()) 
                      : Math.log((speciesAgent.moleculeCount()+1)/phase.volume());        
    }
    
    public double lnProbabilityRatio() {
        double hNew = 0.0;
        if(insert) hNew = -mu + potential.set(phase).calculate(iteratorDirective.set(testMolecule), energy.reset()).sum();
        return -(hNew - hOld)/parentIntegrator.temperature;
    }
    
    public void acceptNotify() {
        if(!insert) testMolecule.sendToReservoir();
    }
    
    public void rejectNotify() {
        if(insert) testMolecule.sendToReservoir();
        else testMolecule.node.setParent((AtomTreeNodeGroup)speciesAgent.node);
    }

    /**
     * Returns an iterator giving molecule that is being added or deleted 
     * in the current or most recent trial.
     */
    public final AtomIterator affectedAtoms(Phase phase) {
        if(this.phase != phase) return AtomIterator.NULL;
        affectedAtomIterator.setBasis(testMolecule);
        affectedAtomIterator.reset();
        return affectedAtomIterator;
    }

    /**
     * Mutator method for the chemical potential of the insertion/deletion species.
     */
    public final void setMu(double mu) {this.mu = mu;}
    /**
     * Accessor method for the chemical potential of th insertion/deletion species.
     */
    public final double getMu() {return mu;}
    /**
     * Indicates that chemical potential has dimensions of energy.
     */
    public final etomica.units.Dimension getMuDimension() {return etomica.units.Dimension.ENERGY;}
/*    
    public static void main(String[] args) {
        etomica.simulations.HsMc2d sim = new etomica.simulations.HsMc2d();
        Simulation.instance = sim;

        MeterNMolecules meterN = new MeterNMolecules();
        etomica.graphics.DisplayBox box = new etomica.graphics.DisplayBox(meterN);
        box.setUpdateInterval(10);
        
        MCMoveInsertDelete mcMoveInsDel = new MCMoveInsertDelete(sim.integrator);
        mcMoveInsDel.setSpecies(sim.species);
        mcMoveInsDel.setMu(-2000.);
        
        sim.integrator(0).setTemperature(1.0);
		                                    
        etomica.graphics.DeviceSlider slider = new etomica.graphics.DeviceSlider(mcMoveInsDel,"mu");
        slider.setMinimum(-20);
        slider.setMaximum(+20);
		Simulation.instance.elementCoordinator.go();

        etomica.graphics.SimulationGraphic.makeAndDisplayFrame(sim);
    }//end of main
// */    
}//end of MCMoveInsertDelete