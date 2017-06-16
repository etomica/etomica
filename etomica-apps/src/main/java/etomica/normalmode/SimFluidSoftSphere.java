/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.action.BoxInflate;
import etomica.action.IAction;
import etomica.action.WriteConfiguration;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.config.ConfigurationFile;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageCollapsing;
import etomica.data.DataPump;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.meter.MeterPressure;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.integrator.mcmove.MCMoveStepTracker;
import etomica.lattice.LatticeCubicFcc;
import etomica.listener.IntegratorListenerAction;
import etomica.potential.P2SoftSphere;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.potential.Potential2SoftSpherical;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.species.SpeciesSpheresMono;

import java.io.File;

/**
 * MC simulation of soft-sphere fluid model  
 * 
 * @author Tai Boon Tan
 */
public class SimFluidSoftSphere extends Simulation {

    private static final long serialVersionUID = 1L;
    private static final String APP_NAME = "Soft Sphere Fluid";
    private static final int PIXEL_SIZE = 10;
    public static String filename;
    public IntegratorMC integrator;
    public SpeciesSpheresMono species;
    public ActivityIntegrate activityIntegrate;
    public Box box;
    public PotentialMaster potentialMaster;
    public SimFluidSoftSphere(Space _space, int numAtoms, double density, double temperature, int exponent) {
        super(_space);

        potentialMaster = new PotentialMaster();

        species = new SpeciesSpheresMono(this, space);
        addSpecies(species);

        box = new Box(space);
        addBox(box);
        box.setNMolecules(species, numAtoms);

        rescaleBox(density);

        integrator = new IntegratorMC(potentialMaster, getRandom(), temperature);
        MCMoveAtom move = new MCMoveAtom(random, potentialMaster, space);
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

        double truncationRadius = box.getBoundary().getBoxSize().getX(0) * 0.495;
        P2SoftSphericalTruncated pTruncated = new P2SoftSphericalTruncated(space, potential, truncationRadius);
        //potentialMaster.lrcMaster().setEnabled(false); //turn off the long-range correction ::updated 7/4/2008

        AtomType sphereType = species.getLeafType();
        potentialMaster.addPotential(pTruncated, new AtomType[]{sphereType, sphereType});

        integrator.setBox(box);
    }

    /**
     * @param args
     */
    public static void main(String[] args) {

    	/*
         *  exponent               freezing density
    	 *  ----------------------------------------
    	 *  	12						0.813
    	 *  	 9						0.943
    	 *  	 6						1.540
    	 *  	 4						3.920
    	 *  Referece: Hoover, Gray & Johnson (JCP 1971) Table IV
    	 */

        // defaults
        int D = 3;
        int nA = 500;
        double density_div_sqrt2 = .215;
        double temperature = 1.0;
        int exponent = 12;
        double freezeDensity = 0.813;
        long simSteps =10000;


        // parse arguments
        if (args.length > 0) {
            density_div_sqrt2 = Double.parseDouble(args[0]);
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
        if (args.length > 5) {
        	freezeDensity = Double.parseDouble(args[5]);
        }

        double density = density_div_sqrt2*Math.sqrt(2);

        int rho60Freeze = (int)Math.round(0.6*freezeDensity*1000);
        double density60Freeze = (double)rho60Freeze/1000;
        int rho90Freeze = (int)Math.round(0.9*freezeDensity*1000);
        double density90Freeze = (double)rho90Freeze/1000;

        System.out.println("Running fluid soft sphere simulation");
        System.out.println(nA + " atoms with exponent " + exponent+" and simulation box density "+density);
        System.out.println("Freezing Density: "+ freezeDensity +" ;Density/sqrt(2): " + density_div_sqrt2);
        System.out.println("isotherm temperature at "+temperature);
        System.out.println(simSteps+ " steps");
        System.out.println("60% Freezing Density: "+ density60Freeze + " ;90% Freezing Density: "+ density90Freeze);


        // construct simulation
        SimFluidSoftSphere sim = new SimFluidSoftSphere(Space.getInstance(D), nA, density, temperature, exponent);

        /*
         * Output the Configuration File Name
         * The format is, e.g. Confign12SS1256
         * 		exponent n = 12
         * 		density/sqrt(2) = 1.256
         */
        if(density_div_sqrt2 < 0.01){
        	filename = "confign"+exponent+"SS000"+(int)Math.round(density_div_sqrt2*1000);

    	} else if (density_div_sqrt2 >= 0.01 && density_div_sqrt2 < 0.1){
        	filename = "confign"+exponent+"SS00"+(int)Math.round(density_div_sqrt2*1000);

        } else if (density_div_sqrt2 >= 0.1 && density_div_sqrt2 < 1.0){
        	filename = "confign"+exponent+"SS0"+(int)Math.round(density_div_sqrt2*1000);

        } else {
        	filename = "confign"+exponent+"SS"+(int)Math.round(density_div_sqrt2*1000);

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
        IntegratorListenerAction pressurePumpListener = new IntegratorListenerAction(pressurePump);
        pressurePumpListener.setInterval(nA/10);
	    sim.integrator.getEventManager().addListener(pressurePumpListener);

        final double d = density;
        final double temp = temperature;
        final int n = exponent;
        IAction pressureCheck = new IAction(){
        	public void actionPerformed(){
        		System.out.println("Var[density/sqrt(2)]: "+ d/(Math.sqrt(2)*Math.pow(temp, 3.0/(double)n))+ " ,Z: "
                        + pressureAverage.getData().getValue(AccumulatorAverage.AVERAGE.index) / (d * temp));

        	}
        };
        IntegratorListenerAction pressureCheckListener = new IntegratorListenerAction(pressureCheck);
        pressureCheckListener.setInterval((int)simSteps/10);
        sim.integrator.getEventManager().addListener(pressureCheckListener);

        /*
        SimulationGraphic simgraphic = new SimulationGraphic(sim, SimulationGraphic.GRAPHIC_ONLY, APP_NAME, sim.space, sim.getController());
        simgraphic.getController().getReinitButton().setPostAction(simgraphic.getPaintAction(sim.box));
        simgraphic.makeAndDisplayFrame(APP_NAME);
        */

        MeterPotentialEnergy meterEnergy = new MeterPotentialEnergy(sim.potentialMaster);
        meterEnergy.setBox(sim.box);

        AccumulatorAverage energyAverage = new AccumulatorAverageCollapsing();
        DataPump energyPump = new DataPump(meterEnergy, energyAverage);
        IntegratorListenerAction energyPumpListener = new IntegratorListenerAction(energyPump);
        energyPumpListener.setInterval(100);
        sim.integrator.getEventManager().addListener(energyPumpListener);


        /*
         * If the density is greater 90% of the freezing density;
         * 	we use the final configuration of the simulation run with 0.9*freezeDensity
         * 	as the initial configuration and scale to the corresponding density.
         *
         */

        File configFileNew = new File(filename + ".pos");

        if(configFileNew.exists()){
    		System.out.println("\n***using "+ configFileNew);
			sim.activityIntegrate.setMaxSteps(10);
        	sim.initializeConfigFromFile(filename);

        } else {

            if (density_div_sqrt2 < density60Freeze){
	        	sim.activityIntegrate.setMaxSteps(simSteps/10);  //simSteps/10

            } else if (density_div_sqrt2 >= density60Freeze && density_div_sqrt2 <= density90Freeze) {
                sim.activityIntegrate.setMaxSteps(simSteps/4);

            } else if (density_div_sqrt2 > density90Freeze){
	        	File configFile09Freeze;
	        	String fileName09Freeze;

                if (density90Freeze < 1.0){
	        		fileName09Freeze = "confign"+exponent+"SS0"+rho90Freeze;
	        		configFile09Freeze = new File(fileName09Freeze+".pos");

                } else {
	        		fileName09Freeze = "confign"+exponent+"SS"+rho90Freeze;
	        		configFile09Freeze = new File(fileName09Freeze+".pos");
	        	}


                if (!configFile09Freeze.exists()){
	        		throw new RuntimeException("Density exceeds 90% of freezing density,\nyou need to have the configuration file from lower density");

                } else {
	        		System.out.println("\n***using "+ configFile09Freeze);
	        		sim.activityIntegrate.setMaxSteps(simSteps/10);
	        		sim.rescaleBox(density90Freeze*Math.sqrt(2));
                    sim.initializeConfigFromFile(fileName09Freeze);
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

        System.out.println("Average Energy: " + energyAverage.getData().getValue(AccumulatorAverage.AVERAGE.index)
                + " ,Error Energy: " + energyAverage.getData().getValue(AccumulatorAverage.ERROR.index));
        System.out.println("Average Pressure: " + pressureAverage.getData().getValue(AccumulatorAverage.AVERAGE.index)
                + " ,Error Pressure: " + pressureAverage.getData().getValue(AccumulatorAverage.ERROR.index));
        System.out.println(" ");

        System.out.println("Var[density/sqrt(2)]: "+ density/(Math.sqrt(2)*Math.pow(temperature, 3/(double)exponent))+ " ,Z: "
                + pressureAverage.getData().getValue(AccumulatorAverage.AVERAGE.index) / (density * temperature));

    }

    public void rescaleBox(double d) {
        BoxInflate boxscaling = new BoxInflate(box, space);
        boxscaling.setTargetDensity(d);
        boxscaling.actionPerformed();
    }

    public void initializeConfigFromFile(String fname) {
        ConfigurationFile config = new ConfigurationFile(fname);
        config.initializeCoordinates(box);
    }

    public void writeConfiguration(String fname) {
        WriteConfiguration writeConfig = new WriteConfiguration(space);
        writeConfig.setBox(box);
        writeConfig.setConfName(fname);
        writeConfig.actionPerformed();
        System.out.println("\n***output configFile: " + fname);
    }
}
