package etomica;

import etomica.units.Dimension;

/**
 * Standard Monte Carlo molecule-displacement trial move.
 *
 * @author David Kofke
 */
 
 /* History of changes
  * 7/9/02 Added energyChange() method
  */
public class MCMoveMolecule extends MCMove {
    
    private final MeterPotentialEnergy energyMeter;
    private final AtomIteratorSinglet affectedAtomIterator = new AtomIteratorSinglet();
    private Atom molecule;
    protected double uOld;
    protected double uNew = Double.NaN;

    public MCMoveMolecule(PotentialMaster potentialMaster) {
        super(potentialMaster);
        energyMeter = new MeterPotentialEnergy(potentialMaster);
        setStepSizeMax(Default.BOX_SIZE);
        setStepSizeMin(0.0);
        setStepSize(0.1*Default.ATOM_SIZE);
        setPerParticleFrequency(true);
        energyMeter.setIncludeLrc(false);
        //set directive to exclude intramolecular contributions to the energy

        //TODO enable meter to do this
        //       iteratorDirective.addCriterion(new IteratorDirective.PotentialCriterion() {
 //           public boolean excludes(Potential p) {return (p instanceof Potential1.Intramolecular);}
 //       });
        setName("MCMoveMolecule");
        
    }
    
    public final Dimension getStepSizeDimension() {return Dimension.LENGTH;}
    public final Dimension getStepSizeMaxDimension() {return Dimension.LENGTH;}
    public final Dimension getStepSizeMinDimension() {return Dimension.LENGTH;}
    

    public boolean doTrial() {
        if(phase.moleculeCount()==0) return false;
        
        molecule = phase.randomMolecule();

        energyMeter.setTarget(molecule);
        uOld = energyMeter.getDataAsScalar(phase);
        molecule.coord.displaceWithin(stepSize);
        uNew = Double.NaN;
        return true;
    }//end of doTrial
    
    public double lnTrialRatio() {return 0.0;}
    
    public double lnProbabilityRatio() {
        energyMeter.setTarget(molecule);
        uNew = energyMeter.getDataAsScalar(phase);
        return -(uNew - uOld)/temperature;
    }
    
    public void acceptNotify() {  /* do nothing */}
    
    public void rejectNotify() {
        molecule.coord.replace();
    }
        
    public final AtomIterator affectedAtoms(Phase phase) {
        if(this.phase != phase) return AtomIterator.NULL;
        affectedAtomIterator.setAtom(molecule);
        affectedAtomIterator.reset();
        return affectedAtomIterator;
    }

    public double energyChange(Phase phase) {return (this.phase == phase) ? uNew - uOld : 0.0;}
    

/*    public static void main(String[] args) {
        
	    IntegratorMC integrator = new IntegratorMC();
	    MCMoveMolecule mcMove = new MCMoveMolecule(integrator);
	    SpeciesSpheres species = new SpeciesSpheres(20,3);
	    Phase phase = new Phase();
	    P2LennardJones potential = new P2LennardJones();
	    Controller controller = new Controller();
	    DisplayPhase display = new DisplayPhase();

		MeterEnergy energy = new MeterEnergy();
		energy.setPhase(phase);
		energy.setHistorying(true);
		energy.setActive(true);		
		energy.getHistory().setNValues(500);		
		DisplayPlot plot = new DisplayPlot();
		plot.setLabel("Energy");
		plot.setDataSources(energy.getHistory());
		
		integrator.setSleepPeriod(2);
		
		Simulation.instance.elementCoordinator.go();
	    
        potential.setIterator(new AtomPairIterator(phase));
        potential.set(species.getAgent(phase));
        
        Simulation.makeAndDisplayFrame(Simulation.instance);
    }//end of main
  */  
}