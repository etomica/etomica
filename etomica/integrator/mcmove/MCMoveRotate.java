package etomica.integrator.mcmove;

import etomica.Atom;
import etomica.AtomIterator;
import etomica.Phase;
import etomica.PotentialMaster;
import etomica.Space;
import etomica.atom.iterator.AtomIteratorSinglet;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.MCMove;
import etomica.space.ICoordinateAngular;
import etomica.space.Orientation;

/**
 * Performs a rotation of an atom (not a molecule) that has an orientation coordinate.
 */

 /* History of changes
  * 7/9/02 Added energyChange() method
  */
  
public class MCMoveRotate extends MCMove {
    
    private MeterPotentialEnergy energyMeter;
    private final AtomIteratorSinglet affectedAtomIterator = new AtomIteratorSinglet();
    private final Orientation oldOrientation;

    private transient Atom molecule;
    private transient double uOld;
    private transient double uNew = Double.NaN;
    private transient Orientation orientation;

    public MCMoveRotate(PotentialMaster potentialMaster, Space space) {
        super(potentialMaster, 1);
        energyMeter = new MeterPotentialEnergy(potentialMaster);
        oldOrientation = space.makeOrientation();
        setStepSizeMax(Math.PI);
        setStepSizeMin(0.0);
        setStepSize(Math.PI/2.0);
        perParticleFrequency = true;
        energyMeter.setIncludeLrc(false);
    }
     
    public boolean doTrial() {
        Phase phase = phases[0];
        if(phase.moleculeCount()==0) {return false;}
        molecule = phase.randomMolecule();

        energyMeter.setTarget(molecule);
        uOld = energyMeter.getDataAsScalar(phase);
        orientation = ((ICoordinateAngular)molecule.coord).orientation(); 
        oldOrientation.E(orientation);  //save old orientation
        orientation.randomRotation(stepSize);
        uNew = Double.NaN;
        return true;
    }//end of doTrial
    
    public double lnTrialRatio() {return 0.0;}
    
    public double lnProbabilityRatio() {
        energyMeter.setTarget(molecule);
        uNew = energyMeter.getDataAsScalar(phases[0]);
        return -(uNew - uOld)/temperature;
    }
    
    public void acceptNotify() {  /* do nothing */}
    
    public void rejectNotify() {
        orientation.E(oldOrientation);
    }

    public double energyChange(Phase phase) {return (this.phases[0] == phase) ? uNew - uOld : 0.0;}
    
    public final AtomIterator affectedAtoms(Phase phase) {
        if(this.phases[0] != phase) return AtomIterator.NULL;
        affectedAtomIterator.setAtom(molecule);
        affectedAtomIterator.reset();
        return affectedAtomIterator;
    }

 /*
    public static void main(String[] args) {
        Default.ATOM_SIZE =1.2;
        Simulation.instance = new Simulation();
        Phase phase1 = new Phase();
        
        IntegratorMC integratorMC1 = new IntegratorMC();
        MCMoveAtom mcmove= new MCMoveAtom(integratorMC1);
        MCMoveRotate mcrotate = new MCMoveRotate(integratorMC1);
        Default.TEMPERATURE = LennardJones.Temperature.UNIT.toSim(1.30);
        //MCMoveVolumeXY mcmovevolume = new MCMoveVolumeXY();
//        MeterPressureByVolumeChange meterp = new MeterPressureByVolumeChange();
	    integratorMC1.setDoSleep(false);

/*	    meterp.setPhase(phase1);
	    meterp.setInflateDimensions(new boolean[] {true, false});
	    meterp.setActive(true);
        etomica.graphics.DisplayPlot plot1 = new etomica.graphics.DisplayPlot();
	    plot1.setDataSources(meterp);
	       
	    plot1.setWhichValue(MeterAbstract.AVERAGE);
* /        integratorMC1.setTemperature(Default.TEMPERATURE);  

	    SpeciesSpheresRotating speciesSpheres1 = new SpeciesSpheresRotating(80);
	  //  SpeciesSpheresRotating speciesSpheres1 = new SpeciesSpheresRotating(200);
	    speciesSpheres1.setDiameter(1.2);
	    
	    P2HardAssociationCone potential = new P2HardAssociationCone();
	    potential.setWellCutoff(1.5*speciesSpheres1.getDiameter());
	    
	    Controller controller1 = new Controller();
	  
	    etomica.graphics.DisplayPhase displayPhase1 = new etomica.graphics.DisplayPhase();
	    MeterEnergy meterEnergy1 = new MeterEnergy();
	    
	    etomica.graphics.DisplayBox box1 = new etomica.graphics.DisplayBox();
	    
	    box1.setMeter(meterEnergy1);
	    box1.setWhichValue(MeterAbstract.AVERAGE);
		
		displayPhase1.setPhase(phase1);
		
		etomica.graphics.DeviceSlider temperatureSlider = new etomica.graphics.DeviceSlider(integratorMC1, "temperature");
	    temperatureSlider.setUnit(new etomica.units.Unit(Kelvin.UNIT));
	    temperatureSlider.setMinimum(50);
	    temperatureSlider.setMaximum(500);
	   
	    phase1.setIntegrator(integratorMC1);
        phase1.boundary().dimensions().setX(0,15);
        
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
        
	//	Simulation.instance.panel().setBackground(java.awt.Color.blue);		
        etomica.graphics.SimulationGraphic.makeAndDisplayFrame(Simulation.instance);
    }//end of main
//   */
}//end of MCMoveRotate
