package etomica;
import java.util.Random;
import etomica.units.*;

public class MCMoveRotate extends MCMove {
    
    private PotentialAgent phasePotential;
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
     
          //need to modify to handle multiple-phase issues
    public void setPhase(Phase p) {
        super.setPhase(p);
        phasePotential = p.potential();
    }

    public void setParentIntegrator(IntegratorMC parent) {
        super.setParentIntegrator(parent);
        oldOrientation = parent.parentSimulation().space().makeOrientation();
    }
     
    public void thisTrial(){   
         
        double uOld, uNew;
       
        if(phase.atomCount==0) {return;}
        int i = (int)(rand.nextDouble()*phase.atomCount);
        Atom a = phase.firstAtom();
        for(int j=i; --j>=0; ) {a = a.nextAtom();}  
            
        phasePotential.calculate(iteratorDirective.set(a), energy.reset());
        uOld = energy.sum();
        Space.Orientation orientation = ((Space.Coordinate.Angular)a.coord).orientation(); 
        oldOrientation.E(orientation);  //save old orientation
        orientation.randomRotation(stepSize);
        phasePotential.calculate(iteratorDirective.set(a), energy.reset());
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
/* 
   public static void main(String[] args) {
        java.awt.Frame f = new java.awt.Frame();   //create a window
        f.setSize(600,350);
        Default.ATOM_SIZE =1.2;
         Simulation.instance = new Simulation(new Space2DCell());
        Phase phase1 = new Phase();
        
        MCMoveAtom mcmove= new MCMoveAtom();
        MCMoveRotate mcrotate = new MCMoveRotate();
        DisplayPlot plot1 = new DisplayPlot();
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
	        plot1.setDataSources(meterp);
	       //plot1.setMeter(meterp);
	       
	        plot1.setWhichValue(MeterAbstract.AVERAGE);
        	integratorMC1.setTemperature(Default.TEMPERATURE);  
        	 System.out.println("temp"+ integratorMC1.getTemperature());
	   // SpeciesSpheresRotating speciesDisks1 = new SpeciesSpheresRotating(80);
	    SpeciesSpheresRotating speciesDisks1 = new SpeciesSpheresRotating(200);
	    speciesDisks1.setDiameter(1.2);
	    // SpeciesDisks    speciesDisk1 = new SpeciesDisks();
	    // speciesDisk1.setNMolecules(200);
	    
	    PotentialAssociationCone potential = new PotentialAssociationCone();
	    potential.setWellCutoff(1.5*speciesDisks1.getDiameter());
	    Potential2 p2 = new P2SimpleWrapper(new PotentialAssociationCone());
//       P2IdealGas p2IdealGas = new P2IdealGas();
	   
	    
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
        
        // Simulation.elementCoordinator = new Simulation.ElementCoordinator.Basic();
		Simulation.instance.elementCoordinator.go(); //invoke this method only after all elements are in place
		                                    //calling it a second time has no effect
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
        f.add(Simulation.instance);         //access the static instance of the simulation to
                                            //display the graphical components
        f.pack();
        f.show();
        phase1.boundary().dimensions().setComponent(0,30);
        f.addWindowListener(new java.awt.event.WindowAdapter() {   //anonymous class to handle window closing
            public void windowClosing(java.awt.event.WindowEvent e) {System.exit(0);}
        });
    }//end of main
  */ 
}
