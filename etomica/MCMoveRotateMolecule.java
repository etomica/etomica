package etomica;
import etomica.units.*;
import etomica.action.AtomActionTransform;

/**
 * Monte Carlo trial that rotates the atoms of a molecule about its first atom.
 *
 * Has a bug, probably associated with incorrect replacement of the molecule when
 * rejecting the trial.
 */
public class MCMoveRotateMolecule extends MCMove {
    
    private final IteratorDirective iteratorDirective = new IteratorDirective(IteratorDirective.BOTH);
    private final AtomIteratorSinglet affectedAtomIterator = new AtomIteratorSinglet();
    
    //not sure that tree iterator will give all leaf atoms of a molecule
    //because it assumes all subgroups have same structure
    private final AtomIterator leafAtomIterator = new AtomIteratorTree();
    
    private transient double uOld;
    private transient Atom molecule;
    private transient Space.Vector r0;
    private transient Space.RotationTensor rotationTensor;

    public MCMoveRotateMolecule(IntegratorMC parentIntegrator) {
        super(parentIntegrator);
        rotationTensor = (Space.RotationTensor)parentIntegrator.parentSimulation().space.makeRotationTensor();
        r0 = parentIntegrator.parentSimulation().space.makeVector();
        setStepSizeMax(Math.PI);
        setStepSizeMin(0.0);
        setStepSize(Math.PI/2.0);
        setPerParticleFrequency(true);
        iteratorDirective.includeLrc = false;
        //set directive to exclude intramolecular contributions to the energy
        iteratorDirective.addCriterion(new IteratorDirective.PotentialCriterion() {
            public boolean excludes(Potential p) {return (p instanceof Potential1.Intramolecular);}
        });
    }
     
    public boolean doTrial() {
        if(phase.moleculeCount()==0) {molecule = null; return false;}
        molecule = phase.randomMolecule();
        leafAtomIterator.setBasis(molecule);

        uOld = potential.set(phase).calculate(iteratorDirective.set(molecule), energy.reset()).sum();
        //update for 3D
        double dTheta = (2*Simulation.random.nextDouble() - 1.0)*stepSize;
        rotationTensor.setAxial(2,dTheta);
       // molecule.coord.transform(molecule.coord.position(), rotationTensor);
        leafAtomIterator.reset();
        r0.E(molecule.node.firstLeafAtom().coord.position());
//        AtomActionTransform.doAction(leafAtomIterator, molecule.coord.position(), rotationTensor);
        AtomActionTransform.doAction(leafAtomIterator, r0, rotationTensor);
        return true;
    }//end of doTrial
    
    public double lnTrialRatio() {return 0.0;}
    
    public double lnProbabilityRatio() {
        double uNew = potential.calculate(iteratorDirective, energy.reset()).sum();
        return -(uNew - uOld)/parentIntegrator.temperature;
    }
    
    public void acceptNotify() {  /* do nothing */}
    
    public void rejectNotify() {
        rotationTensor.invert();
        leafAtomIterator.reset();
        AtomActionTransform.doAction(leafAtomIterator, r0, rotationTensor);
    }
 
    public final AtomIterator affectedAtoms(Phase phase) {
        if(this.phase != phase) return AtomIterator.NULL;
        affectedAtomIterator.setBasis(molecule);
        affectedAtomIterator.reset();
        return affectedAtomIterator;
    }

/*   public static void main(String[] args) {
        Default.ATOM_SIZE =1.2;
        Simulation.instance = new Simulation(/*new Space2DCell()* /);
        Phase phase1 = new Phase();
        
        IntegratorMC integratorMC1 = new IntegratorMC();
        MCMoveAtom mcmove= new MCMoveAtom(integratorMC1);
        MCMoveRotate mcrotate = new MCMoveRotate(integratorMC1);
         Default.TEMPERATURE = LennardJones.Temperature.UNIT.toSim(1.30);
        //MCMoveVolumeXY mcmovevolume = new MCMoveVolumeXY();
        MeterPressureByVolumeChange meterp = new MeterPressureByVolumeChange();
	    integratorMC1.setDoSleep(false);

	    meterp.setPhase(phase1);
	    meterp.setInflateDimensions(new boolean[] {true, false});
	    meterp.setActive(true);
        DisplayPlot plot1 = new DisplayPlot();
	    plot1.setDataSources(meterp);
	       
	    plot1.setWhichValue(MeterAbstract.AVERAGE);
        integratorMC1.setTemperature(Default.TEMPERATURE);  

	    SpeciesSpheresRotating speciesSpheres1 = new SpeciesSpheresRotating(80);
	  //  SpeciesSpheresRotating speciesSpheres1 = new SpeciesSpheresRotating(200);
	    speciesSpheres1.setDiameter(1.2);
	    
	    P2HardAssociationCone potential = new P2HardAssociationCone();
	    potential.setWellCutoff(1.5*speciesSpheres1.getDiameter());
	    
	    Controller controller1 = new Controller();
	  
	    DisplayPhase displayPhase1 = new DisplayPhase();
	    MeterEnergy meterEnergy1 = new MeterEnergy();
	    
	    DisplayBox box1 = new DisplayBox();
	    
	    box1.setMeter(meterEnergy1);
	    box1.setWhichValue(MeterAbstract.AVERAGE);
		
		displayPhase1.setPhase(phase1);
		
		DeviceSlider temperatureSlider = new DeviceSlider(integratorMC1, "temperature");
	    temperatureSlider.setUnit(new etomica.units.Unit(Kelvin.UNIT));
	    temperatureSlider.setMinimum(50);
	    temperatureSlider.setMaximum(500);
	   
	    phase1.setIntegrator(integratorMC1);
        phase1.boundary().dimensions().setComponent(0,15);
        
		Simulation.instance.elementCoordinator.go(); 
        //phase1.integrator().setTemperature();
        //meterEnergy1.setPhase(phase1);
        //meterpxx.setPhase(phase1);
        //phase1.setIntegrator(integratorMC1);
       // integratorMC1.addIntervalListener(box1);
       // integratorMC1.addIntervalListener(meterp);
        //integratorMC1.setTemperature(Kelvin.UNIT.toSim(10));
        //integratorMC1.addIntervalListener(box2);
//        integratorMC1.addIntervalListener(phase1);
     //  controller1.add(integratorMC1);
        
		Simulation.instance.panel().setBackground(java.awt.Color.blue);		
	//	((simulate.Space2DCell.CellListIteratorFactory)phase1.iteratorFactory()).setNeighborDistance(1.2*Default.ATOM_SIZE);
    //    ((simulate.Space2DCell.CellListIteratorFactory)phase1.iteratorFactory()).setNCells(6,10);
        Simulation.makeAndDisplayFrame(Simulation.instance);
    }//end of main
   */
}
