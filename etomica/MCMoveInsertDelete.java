package etomica;

/**
 * Elementary Monte Carlo move in which a molecule of a specified species is
 * inserted into or removed from a phase.
 *
 * @author David Kofke
 */
 
 /* History of changes
  * 07/09/02 (DAK) Added energyChange() method
  * 09/19/02 (DAK) Minor change in doTrial for case were deleting with N = 0
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
    private double uOld;
    private double uNew = Double.NaN;
    private boolean insert;
    private AtomReservoir reservoir;

    public MCMoveInsertDelete(IntegratorMC parent) {
        super(parent);
        setStepSizeMax(1.0);
        setStepSizeMin(0.0);
        setStepSize(0.10);
        setMu(0.0);
        setTunable(false);
        iteratorDirective.includeLrc = true;
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
            uOld = 0.0;
            testMolecule = species.moleculeFactory().makeAtom((AtomTreeNodeGroup)speciesAgent.node);
            testMolecule.coord.translateTo(phase.randomPosition());
        } else {//delete
            if(speciesAgent.moleculeCount() == 0) {
                testMolecule = null;//added this line 09/19/02
                return false;
            }
            testMolecule = speciesAgent.randomMolecule();
            uOld = potential.calculate(phase, iteratorDirective.set(testMolecule), energy.reset()).sum();
            reservoir.addAtom(testMolecule);
        } 
        uNew = Double.NaN;
        return true;
    }//end of doTrial
    
    public double lnTrialRatio() {//note that moleculeCount() gives the number of molecules after the trial is attempted
        return insert ? Math.log(phase.volume()/speciesAgent.moleculeCount()) 
                      : Math.log((speciesAgent.moleculeCount()+1)/phase.volume());        
    }
    
    public double lnProbabilityRatio() {
        if(insert) {
            uNew = potential.calculate(phase, iteratorDirective.set(testMolecule), energy.reset()).sum();
            return (+mu - uNew)/parentIntegrator.temperature;
        } else {//delete
            uNew = 0.0;
            return (-mu + uOld)/parentIntegrator.temperature;
        }
    }
    
    public void acceptNotify() {
        if(!insert) testMolecule.sendToReservoir();
    }
    
    public void rejectNotify() {
        if(insert) testMolecule.sendToReservoir();
        else testMolecule.node.setParent((AtomTreeNodeGroup)speciesAgent.node);
    }
    
    public double energyChange(Phase phase) {return (this.phase == phase) ? uNew - uOld : 0.0;}

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
///*    
    public static void main(String[] args) {
        Default.TRUNCATE_POTENTIALS = false;
//        etomica.simulations.HsMc2d sim = new etomica.simulations.HsMc2d();
//		MeterNMolecules meterN = new MeterNMolecules();
//		etomica.graphics.DisplayBox box = new etomica.graphics.DisplayBox((DatumSource)meterN);
//		box.setUpdateInterval(10);
//        
//		MCMoveInsertDelete mcMoveInsDel = new MCMoveInsertDelete(sim.integrator);
//		mcMoveInsDel.setSpecies(sim.species);
//		mcMoveInsDel.setMu(-2000.);
//        
//		sim.integrator(0).setTemperature(1.0);
//		                                    
//		etomica.graphics.DeviceSlider slider = new etomica.graphics.DeviceSlider(mcMoveInsDel,"mu");
//		slider.setMinimum(-20);
//		slider.setMaximum(+20);
//		Simulation.instance.elementCoordinator.go();
		Simulation sim = new etomica.graphics.SimulationGraphic();
		Controller controller = new Controller();

		SpeciesSpheresMono species = new SpeciesSpheresMono();
		
		Phase phase1 = new Phase();
		Phase phase2 = new Phase();
		IntegratorMC integrator1 = new IntegratorMC();
		IntegratorMC integrator2 = new IntegratorMC();
		phase1.setIntegrator(integrator1);
		phase2.setIntegrator(integrator2);
		integrator1.setTemperature(1.0);
		integrator2.setTemperature(1.0);
		
		P2HardSphere potential = new P2HardSphere();
		potential.setSpecies(species);
		
		sim.elementCoordinator.go();
		
		MCMoveAtom moveAtom1 = new MCMoveAtom(integrator1);
		MCMoveAtom moveAtom2 = new MCMoveAtom(integrator2);
		MCMoveInsertDelete mcMoveInsDel1 = new MCMoveInsertDelete(integrator1);
		MCMoveInsertDelete mcMoveInsDel2 = new MCMoveInsertDelete(integrator2);
		mcMoveInsDel1.setSpecies(species);
		mcMoveInsDel2.setSpecies(species);		
		mcMoveInsDel1.setMu(-2000.);
		mcMoveInsDel2.setMu(0.);
		
		MeterNMolecules meterN1 = new MeterNMolecules();
		meterN1.setPhase(phase1);
		etomica.graphics.DisplayBox box1 = new etomica.graphics.DisplayBox((DatumSource)meterN1);
		box1.setUpdateInterval(10);
		MeterNMolecules meterN2 = new MeterNMolecules();
		meterN2.setPhase(phase2);
		etomica.graphics.DisplayBox box2 = new etomica.graphics.DisplayBox((DatumSource)meterN2);
		box1.setUpdateInterval(10);
        
		etomica.graphics.DeviceSlider slider = new etomica.graphics.DeviceSlider(mcMoveInsDel1,"mu");
		slider.setMinimum(-20);
		slider.setMaximum(+20);
		
		etomica.graphics.DisplayPhase display1 = new etomica.graphics.DisplayPhase();
		etomica.graphics.DisplayPhase display2 = new etomica.graphics.DisplayPhase();
		display1.setPhase(phase1);
		display2.setPhase(phase2);
		Simulation.instance.elementCoordinator.go();

        etomica.graphics.SimulationGraphic.makeAndDisplayFrame(sim);
    }//end of main
// */   
}//end of MCMoveInsertDelete