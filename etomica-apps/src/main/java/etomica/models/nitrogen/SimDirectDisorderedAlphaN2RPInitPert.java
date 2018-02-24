/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.nitrogen;

import etomica.action.WriteConfiguration;
import etomica.action.activity.ActivityIntegrate;
import etomica.box.Box;
import etomica.config.ConfigurationFile;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.DataPump;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.IntegratorMC;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisCubicFcc;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.integrator.IntegratorListenerAction;
import etomica.nbr.list.molecule.PotentialMasterListMolecular;
import etomica.normalmode.BasisBigCell;
import etomica.normalmode.MCMoveMoleculeCoupled;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.species.ISpecies;
import etomica.species.Species;
import etomica.units.Kelvin;

import java.io.File;

/**
 * Direct Sampling for Rotational Perturbation
 * from NO rotational energy to constraint-angle rotational energy
 * 
 * @author Tai Boon Tan
 */
public class SimDirectDisorderedAlphaN2RPInitPert extends Simulation {

    public SimDirectDisorderedAlphaN2RPInitPert(Space space, int numMolecules, double density, double temperature, double angle, long numSteps) {
        super(space);

        int nC = (int) Math.pow(numMolecules / 3.999999999, 1.0 / 3.0);
        double a = Math.pow(4.0 / density, 1.0 / 3.0);
        System.out.println("Unit Cell Length, a: " + a);

        Basis basisFCC = new BasisCubicFcc();
        Basis basis = new BasisBigCell(space, basisFCC, new int[]{nC, nC, nC});

        Species species = new SpeciesN2(space);
        addSpecies(species);

        Boundary boundary = new BoundaryRectangularPeriodic(space, nC * a);
        box = this.makeBox(boundary);
        box.setNMolecules(species, numMolecules);

        int[] nCells = new int[]{1, 1, 1};
        Primitive primitive = new PrimitiveCubic(space, nC * a);

        CoordinateDefinitionNitrogen coordinateDef = new CoordinateDefinitionNitrogen(this, box, primitive, basis, space);
        coordinateDef.setIsAlpha();
        coordinateDef.setOrientationVectorAlpha(space);
        coordinateDef.initializeCoordinates(nCells);

        double rcScale = 0.475;
        double rC = box.getBoundary().getBoxSize().getX(0) * rcScale;
        System.out.println("Truncation Radius (" + rcScale + " Box Length): " + rC);
        P2Nitrogen potential = new P2Nitrogen(space, rC);
        potential.setBox(box);

        PotentialMasterListMolecular potentialMaster = new PotentialMasterListMolecular(this, space);
        potentialMaster.addPotential(potential, new ISpecies[]{species, species});

        MCMoveMoleculeCoupled move = new MCMoveMoleculeCoupled(potentialMaster, getRandom(), space);
        move.setBox(box);
        move.setPotential(potential);
        move.setDoExcludeNonNeighbors(true);

        IntegratorMC integrator = new IntegratorMC(potentialMaster, getRandom(), temperature, box);
        integrator.getMoveManager().addMCMove(move);

        int cellRange = 6;
        potentialMaster.setRange(rC);
        potentialMaster.setCellRange(cellRange);
        potentialMaster.getNeighborManager(box).reset();
        potential.setRange(Double.POSITIVE_INFINITY);

        MeterPotentialEnergy meterPE = new MeterPotentialEnergy(potentialMaster, box);

        double latticeEnergy = meterPE.getDataAsScalar();
        System.out.println("lattice energy (sim unit): " + latticeEnergy);
        System.out.println("lattice energy per mol(sim unit): " + latticeEnergy / numMolecules);

        MeterDirectInitPert meterDirectInitPert = new MeterDirectInitPert(integrator, meterPE, coordinateDef, getRandom());
        meterDirectInitPert.setConstraintAngle(angle);

        int numBlock = 100;
        int interval = numMolecules;
        long blockSize = numSteps / (numBlock * interval);
        if (blockSize == 0) blockSize = 1;
        System.out.println("block size " + blockSize + " interval " + interval);
        boltzmannAverage = new AccumulatorAverageFixed(blockSize);

        DataPump boltzmannPump = new DataPump(meterDirectInitPert, boltzmannAverage);
        IntegratorListenerAction boltzmannPumpListener = new IntegratorListenerAction(boltzmannPump, 100);
        integrator.getEventManager().addListener(boltzmannPumpListener);

        activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);
    }

    public void initializeConfigFromFile(String fname){
        ConfigurationFile config = new ConfigurationFile(fname);
        config.initializeCoordinates(box);
    }
    
    public void writeConfiguration(String fname){
        WriteConfiguration writeConfig = new WriteConfiguration(space);
        writeConfig.setBox(box);
        writeConfig.setDoApplyPBC(false);
        writeConfig.setConfName(fname);
        writeConfig.actionPerformed();
        System.out.println("\n***output configFile: "+ fname);
    }
    
    
    /**
     * @param args filename containing simulation parameters
     * @see SimDirectBetaN2RPAngleToNoAngle.SimOverlapParam
     */
    public static void main(String[] args) {

        double temperature = 50; //in UNIT KELVIN
        double density = 0.0220;
        double angle =1;
        long numSteps = 100000;
        int numMolecules = 864;

    	if(args.length > 0){
			temperature = Double.parseDouble(args[0]);
		}
		if(args.length > 1){
			density = Double.parseDouble(args[1]);
		}
		if(args.length > 2){
			angle = Double.parseDouble(args[2]);
		}
		if(args.length > 3){
			numSteps = Long.parseLong(args[3]);
		}
		if(args.length > 4){
			numMolecules = Integer.parseInt(args[4]);
		}
		String configFileName = "configT"+temperature+"_d"+density;
	        
        System.out.println("Running disordered alpha-phase Nitrogen RP direct sampling simulation");
        System.out.println(numMolecules+" molecules at density "+density+" and temperature "+temperature + " K");
        System.out.println("perturbing from no rotational energy into rotational dof of angle= "+angle);
        System.out.println("with numStep of "+ numSteps);
        
        SimDirectDisorderedAlphaN2RPInitPert sim = new SimDirectDisorderedAlphaN2RPInitPert(Space.getInstance(3), numMolecules, density, Kelvin.UNIT.toSim(temperature), angle, numSteps);

        
        File configFile = new File(configFileName+".pos");
        if(configFile.exists()){
			System.out.println("\n***Initialize coordinate from "+ configFile);
        	sim.initializeConfigFromFile(configFileName);
		} else {
			long equiStep = (numMolecules*numSteps/1000);
	        System.out.println("\nEquilibration step: " + equiStep);
	        sim.activityIntegrate.setMaxSteps(equiStep);
	        sim.getController().actionPerformed();     
	        System.out.println("Equilibration finished");
	        sim.getController().reset();
		}
        
        long startTime = System.currentTimeMillis();
        System.out.println("Start Time: " + startTime);
       
        sim.activityIntegrate.setMaxSteps(numSteps);
        sim.getController().actionPerformed();

        sim.writeConfiguration(configFileName);
        double average = sim.boltzmannAverage.getData().getValue(sim.boltzmannAverage.AVERAGE.index);
        double error = sim.boltzmannAverage.getData().getValue(sim.boltzmannAverage.ERROR.index);
        double blockCorrelation = sim.boltzmannAverage.getData().getValue(sim.boltzmannAverage.BLOCK_CORRELATION.index);
        
        System.out.println("blockCorrelation: " + blockCorrelation);
        System.out.println("boltzmann average: " + average + " ;err: " + error +" ;errC: "+ error*Math.sqrt((1+blockCorrelation)/(1-blockCorrelation)));
        
        long endTime = System.currentTimeMillis();
        System.out.println("\nEnd Time: " + endTime);
        System.out.println("Time taken (s): " + (endTime - startTime)/1000);
       
    }

    private static final long serialVersionUID = 1L;
    protected ActivityIntegrate activityIntegrate;
    protected AccumulatorAverageFixed boltzmannAverage;
    protected Box box;

}
