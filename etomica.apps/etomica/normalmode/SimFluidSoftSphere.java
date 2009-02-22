package etomica.normalmode;

import java.io.File;

import etomica.action.BoxInflate;
import etomica.action.WriteConfiguration;
import etomica.action.activity.ActivityIntegrate;
import etomica.api.IAction;
import etomica.api.IAtomType;
import etomica.api.IBox;
import etomica.box.Box;
import etomica.config.ConfigurationFile;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageCollapsing;
import etomica.data.DataPump;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.meter.MeterPressure;
import etomica.data.types.DataGroup;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.integrator.mcmove.MCMoveStepTracker;
import etomica.lattice.LatticeCubicFcc;
import etomica.potential.P2SoftSphere;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.potential.Potential2SoftSpherical;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.species.SpeciesSpheresMono;

/**
 * MC simulation of soft-sphere fluid model  
 * 
 * @author Tai Boon Tan
 */
public class SimFluidSoftSphere extends Simulation {

    public SimFluidSoftSphere(Space _space, int numAtoms, double density, double temperature, int exponent) {
        super(_space, true);

        potentialMaster = new PotentialMaster();

        species = new SpeciesSpheresMono(this, space);
        getSpeciesManager().addSpecies(species);

        box = new Box(space);
        addBox(box);
        box.setNMolecules(species, numAtoms);
   
        rescaleBox(density);
        
        integrator = new IntegratorMC(potentialMaster, getRandom(), temperature);
        MCMoveAtom move = new MCMoveAtom(this, potentialMaster, space);
        move.setStepSize(0.2);
       // move.setStepSizeMax(0.5);
        integrator.getMoveManager().addMCMove(move);
        ((MCMoveStepTracker)move.getTracker()).setNoisyAdjustment(true);
        ((MCMoveStepTracker)move.getTracker()).setAdjustInterval(10);
        
        activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);

       	ConfigurationLattice config = new ConfigurationLattice(new LatticeCubicFcc(space), space);
       	config.initializeCoordinates(box);
       	
        Potential2SoftSpherical potential = new P2SoftSphere(space);
        
        double truncationRadius = box.getBoundary().getDimensions().x(0) * 0.495;
        P2SoftSphericalTruncated pTruncated = new P2SoftSphericalTruncated(space, potential, truncationRadius);
        //potentialMaster.lrcMaster().setEnabled(false); //turn off the long-range correction ::updated 7/4/2008 
        
        IAtomType sphereType = species.getLeafType();
        potentialMaster.addPotential(pTruncated, new IAtomType[] {sphereType, sphereType});
      
        integrator.setBox(box);
    }

    public void rescaleBox(double d){
        BoxInflate boxscaling = new BoxInflate(box, space);
        boxscaling.setTargetDensity(d);
        boxscaling.actionPerformed();
    }
    
    public void initializeConfigFromFile(String fname){
        ConfigurationFile config = new ConfigurationFile(fname);
        config.initializeCoordinates(box);
    }
    
    public void writeConfiguration(String fname){
        WriteConfiguration writeConfig = new WriteConfiguration(space);
        writeConfig.setBox(box);
        writeConfig.setConfName(fname);
        writeConfig.actionPerformed();
        System.out.println("***output configFile: "+ fname);
    }
    
    /**
     * @param args
     */
    public static void main(String[] args) {

        // defaults
        int D = 3;
        int nA = 500;
        double density = 0.05;
        double temperature = 1.0;
        int exponent = 12 ;
        long simSteps =10000;
        
        
        // parse arguments
        if (args.length > 0) {
            density = Double.parseDouble(args[0]);
        }
        if (args.length > 1) {
            simSteps = Long.parseLong(args[1]);
        }
        if (args.length > 2) {
            nA = Integer.parseInt(args[2]);
        }
        if (args.length > 3) {
            temperature = Double.parseDouble(args[3]);
        }
        if (args.length > 4) {
        	exponent = Integer.parseInt(args[4]);
        }
        
        System.out.println("Running fluid soft sphere simulation");
        System.out.println(nA + " atoms with exponent " + exponent+" and density "+density);
        System.out.println("isotherm temperature at "+temperature);
        System.out.println(simSteps+ " steps");

        // construct simulation
        SimFluidSoftSphere sim = new SimFluidSoftSphere(Space.getInstance(D), nA, density, temperature, exponent);
        
        if (density < 0.1){
        	filename = "ConfigSS00"+(int)Math.round(density*100);
        	
        } else if (density >= 0.1 && density < 1.0){
        	filename = "ConfigSS0"+(int)Math.round(density*100);
        	
        } else {
        	filename = "ConfigSS"+(int)Math.round(density*100);
        	
        }
        
        /*
         * Simulation Graphic
         */
        /*
        SimulationGraphic simGraphic = new SimulationGraphic(sim, APP_NAME, 1, sim.space,sim.getController());
        Pixel pixel = new Pixel(10);
        simGraphic.getDisplayBox(sim.box).setPixelUnit(pixel);
        
        simGraphic.getController().getReinitButton().setPostAction(simGraphic.getPaintAction(sim.box));
        simGraphic.getDisplayBox(sim.box).setPixelUnit(new Pixel(PIXEL_SIZE));
        simGraphic.makeAndDisplayFrame(APP_NAME);
        simGraphic.getDisplayBox(sim.box).repaint();
        */
        
       // MeterWidomInsertion meterInsertion = new MeterWidomInsertion(Space.getInstance(D), sim.getRandom());
        //meterInsertion.setIntegrator(sim.integrator);
        //meterInsertion.setSpecies(sim.species);
        //meterInsertion.setNInsert();
       /* 
        AccumulatorAverage insertionAverage = new AccumulatorAverageCollapsing();
        DataPump insertionPump = new DataPump(meterInsertion, insertionAverage);
        sim.integrator.addIntervalAction(insertionPump);
        sim.integrator.setActionInterval(insertionPump, 1);
        */


        MeterPressure meterPressure = new MeterPressure(sim.space);
        meterPressure.setIntegrator(sim.integrator);
        
        final AccumulatorAverage pressureAverage = new AccumulatorAverageCollapsing();
	    DataPump pressurePump = new DataPump(meterPressure, pressureAverage);
	    sim.integrator.addIntervalAction(pressurePump);
	    sim.integrator.setActionInterval(pressurePump, nA/10);
        
        final double d = density; 
        final double temp = temperature;
        IAction pressureCheck = new IAction(){
        	public void actionPerformed(){
        		System.out.println("Var[density/sqrt(2)]: "+ d/(Math.sqrt(2)*Math.pow(temp, 0.25))+ " ,Z: "
                		+((DataGroup)pressureAverage.getData()).getValue(AccumulatorAverage.StatType.AVERAGE.index)/(d*temp));
                
        	}
        };
        
		
        sim.integrator.addIntervalAction(pressureCheck);
        sim.integrator.setActionInterval(pressureCheck, (int)simSteps/10);
        
        /*
        SimulationGraphic simgraphic = new SimulationGraphic(sim, SimulationGraphic.GRAPHIC_ONLY, APP_NAME, sim.space, sim.getController());
        simgraphic.getController().getReinitButton().setPostAction(simgraphic.getPaintAction(sim.box));
        simgraphic.makeAndDisplayFrame(APP_NAME);
        */
                
        MeterPotentialEnergy meterEnergy = new MeterPotentialEnergy(sim.potentialMaster);
        meterEnergy.setBox(sim.box);
              
        AccumulatorAverage energyAverage = new AccumulatorAverageCollapsing();
        DataPump energyPump = new DataPump(meterEnergy, energyAverage);
        sim.integrator.addIntervalAction(energyPump);
        sim.integrator.setActionInterval(energyPump, 100);
        
        
        /*
         * Note: 	density= 1.05 @ kT = 1.0
         * 			var = 0.74246212
         */
        
        File configFileNew = new File(filename+".pos"); 
       
        if(configFileNew.exists()){
    		System.out.println("\n***using "+ configFileNew);
			sim.activityIntegrate.setMaxSteps(10);
        	sim.initializeConfigFromFile(filename);
        
        } else {
        	        
	        if (density < 0.85){
	        	sim.activityIntegrate.setMaxSteps(simSteps/10);  //simSteps/10
	        
	        } else if(density >= 0.85 && density <= 1.05){	
	        	sim.activityIntegrate.setMaxSteps(simSteps/4);
	        	
	        } else if (density > 1.05){
	        	File configFile105 = new File("ConfigSS105.pos_new");
	        	
	        	if (!configFile105.exists()){
	        		throw new RuntimeException("you need to have the configuration file from lower density");
	        	
	        	} else {
	        		System.out.println("\n***using "+ configFile105);
	        		sim.activityIntegrate.setMaxSteps(simSteps/10);
	        		sim.rescaleBox(1.05);
	        		sim.initializeConfigFromFile("ConfigSS105");       	
	        		sim.rescaleBox(density);
	        	
	        	}
	        
	        }
        }
        
        sim.getController().actionPerformed();
        System.out.println("equilibrated");
        
        sim.integrator.getMoveManager().setEquilibrating(false);
        pressureAverage.reset();
        sim.getController().reset();
        
        sim.activityIntegrate.setMaxSteps(simSteps);
        sim.getController().actionPerformed();
        /*
        double insertionScalar = ((DataGroup)insertionAverage.getData()).getValue(AccumulatorAverage.StatType.AVERAGE.index);
        System.out.println("Average insertion scalar: "+ insertionScalar);
        System.out.println("Error insertion scalar: "+ ((DataGroup)insertionAverage.getData()).getValue(AccumulatorAverage.StatType.ERROR.index));
        System.out.println("Liquid Gibbs free energy: " + -temperature*Math.log(insertionScalar));
        System.out.println(" ");
        */
        
        sim.writeConfiguration(filename);
        
        System.out.println("Average Energy: "+ ((DataGroup)energyAverage.getData()).getValue(AccumulatorAverage.StatType.AVERAGE.index)
        					+ " ,Error Energy: "+ ((DataGroup)energyAverage.getData()).getValue(AccumulatorAverage.StatType.ERROR.index));
        System.out.println("Average Pressure: "+ ((DataGroup)pressureAverage.getData()).getValue(AccumulatorAverage.StatType.AVERAGE.index)
        					+ " ,Error Pressure: "+ ((DataGroup)pressureAverage.getData()).getValue(AccumulatorAverage.StatType.ERROR.index));
        System.out.println(" ");
        
        System.out.println("Var[density/sqrt(2)]: "+ density/(Math.sqrt(2)*Math.pow(temperature, 0.25))+ " ,Z: "
        		+((DataGroup)pressureAverage.getData()).getValue(AccumulatorAverage.StatType.AVERAGE.index)/(density*temperature));

    }

    private static final long serialVersionUID = 1L;
    private static final String APP_NAME = "Soft Sphere Fluid";
    private static final int PIXEL_SIZE = 10;
    public IntegratorMC integrator;
    public SpeciesSpheresMono species;
    public ActivityIntegrate activityIntegrate;
    public IBox box;
    public PotentialMaster potentialMaster;
    public static String filename;
}