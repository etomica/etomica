package etomica;
import java.util.Random;
import etomica.units.*;

public class MCMoveRotate extends MCMove {
    
    private final IteratorDirective iteratorDirective = new IteratorDirective(IteratorDirective.BOTH);
    private final PotentialCalculation.EnergySum energy = new PotentialCalculation.EnergySum();
    private final Random rand = new Random();
    private Space.Orientation oldOrientation;

    public MCMoveRotate() {
        setStepSizeMax(Math.PI);
        setStepSizeMin(0.0);
        setStepSize(Math.PI/2.0);
        setPerParticleFrequency(true);
    }

    public void setParentIntegrator(IntegratorMC parent) {
        super.setParentIntegrator(parent);
        oldOrientation = parent.parentSimulation().space().makeOrientation();
    }
     
    public void thisTrial(){   
         
        double uOld, uNew;
       
        if(phase.atomCount()==0) {return;}
        int i = (int)(rand.nextDouble()*phase.atomCount());
        Atom a = phase.firstAtom();
        for(int j=i; --j>=0; ) {a = a.nextAtom();}  
            
        potential.set(phase).calculate(iteratorDirective.set(a), energy.reset());
        uOld = energy.sum();
        Space.Orientation orientation = ((Space.Coordinate.Angular)a.coord).orientation(); 
        oldOrientation.E(orientation);  //save old orientation
        orientation.randomRotation(stepSize);
        potential.set(phase).calculate(iteratorDirective.set(a), energy.reset());
        uNew = energy.sum();
        if(uNew < uOld) {   //accept
            nAccept++;
            return;
        }
        if(uNew >= Double.MAX_VALUE ||  //reject
            Math.exp(-(uNew-uOld)/parentIntegrator.temperature) < rand.nextDouble()) {
                // restore the old values of orientation.
            orientation.E(oldOrientation);
            return;
        }
        nAccept++;   //accept
    }
 
   public static void main(String[] args) {
        Default.ATOM_SIZE =1.2;
        Simulation.instance = new Simulation(/*new Space2DCell()*/);
        Phase phase1 = new Phase();
        
        MCMoveAtom mcmove= new MCMoveAtom();
        MCMoveRotate mcrotate = new MCMoveRotate();
         Default.TEMPERATURE = LennardJones.Temperature.UNIT.toSim(1.30);
        //MCMoveVolumeXY mcmovevolume = new MCMoveVolumeXY();
        MeterPressureByVolumeChange meterp = new MeterPressureByVolumeChange();
        IntegratorMC integratorMC1 = new IntegratorMC();
	    integratorMC1.add(mcmove);
	    integratorMC1.setDoSleep(false);

	    integratorMC1.add(mcrotate);
	     //integratorMC1.add(mcmovevolume);
	    meterp.setPhase(phase1);
	    meterp.setInflateDimension(0);
	    meterp.setActive(true);
        DisplayPlot plot1 = new DisplayPlot();
	    plot1.setDataSources(meterp);
	       
	    plot1.setWhichValue(MeterAbstract.AVERAGE);
        integratorMC1.setTemperature(Default.TEMPERATURE);  

	    SpeciesSpheresRotating speciesDisks1 = new SpeciesSpheresRotating(80);
	  //  SpeciesSpheresRotating speciesDisks1 = new SpeciesSpheresRotating(200);
	    speciesDisks1.setDiameter(1.2);
	    
	    P2HardAssociationCone potential = new P2HardAssociationCone();
	    potential.setWellCutoff(1.5*speciesDisks1.getDiameter());
	    
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
        
		Simulation.instance.setBackground(java.awt.Color.blue);		
	//	((simulate.Space2DCell.CellListIteratorFactory)phase1.iteratorFactory()).setNeighborDistance(1.2*Default.ATOM_SIZE);
    //    ((simulate.Space2DCell.CellListIteratorFactory)phase1.iteratorFactory()).setNCells(6,10);
        Simulation.makeAndDisplayFrame(Simulation.instance);
    }//end of main
   
}
