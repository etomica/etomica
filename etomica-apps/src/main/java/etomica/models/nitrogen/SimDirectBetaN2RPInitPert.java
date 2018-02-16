/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.nitrogen;

import etomica.action.WriteConfiguration;
import etomica.action.activity.ActivityIntegrate;
import etomica.box.Box;
import etomica.box.BoxAgentManager;
import etomica.config.ConfigurationFile;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.DataPump;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.IntegratorMC;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisHcp;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveHexagonal;
import etomica.integrator.IntegratorListenerAction;
import etomica.nbr.cell.molecule.NeighborCellManagerMolecular;
import etomica.nbr.list.molecule.BoxAgentSourceCellManagerListMolecular;
import etomica.nbr.list.molecule.NeighborListManagerSlantyMolecular;
import etomica.nbr.list.molecule.PotentialMasterListMolecular;
import etomica.normalmode.BasisBigCell;
import etomica.normalmode.MCMoveMoleculeCoupled;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryDeformablePeriodic;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.ISpecies;
import etomica.units.Degree;
import etomica.units.Kelvin;

import java.io.File;

/**
 * Direct Sampling for Rotational Perturbation
 * from NO rotational energy to constraint-angle rotational energy
 * 
 * @author Tai Boon Tan
 */
public class SimDirectBetaN2RPInitPert extends Simulation {

    public SimDirectBetaN2RPInitPert(Space space, int numMolecules, double density, double temperature, double angle, long numSteps) {
		super(space);

		BoxAgentSourceCellManagerListMolecular boxAgentSource = new BoxAgentSourceCellManagerListMolecular(this, null, space);
		BoxAgentManager<NeighborCellManagerMolecular> boxAgentManager = new BoxAgentManager<NeighborCellManagerMolecular>(boxAgentSource, this);

		SpeciesN2 species = new SpeciesN2(space);
		addSpecies(species);

		box = new Box(space);
		addBox(box);
		box.setNMolecules(species, numMolecules);

		double ratio = 1.631;
		double aDim = Math.pow(4.0 / (Math.sqrt(3.0) * ratio * density), 1.0 / 3.0);
		double cDim = aDim * ratio;
		System.out.println("\naDim: " + aDim + " ;cDim: " + cDim);
		int nC = (int) Math.pow(numMolecules / 1.999999999, 1.0 / 3.0);

		Basis basisHCP = new BasisHcp();
		BasisBigCell basis = new BasisBigCell(space, basisHCP, new int[]{nC, nC, nC});

		Vector[] boxDim = new Vector[3];
		boxDim[0] = Vector.of(nC * aDim, 0, 0);
		boxDim[1] = Vector.of(
				-nC * aDim * Math.cos(Degree.UNIT.toSim(60)),
				nC * aDim * Math.sin(Degree.UNIT.toSim(60)),
				0
		);
		boxDim[2] = Vector.of(0, 0, nC * cDim);

		int[] nCells = {1, 1, 1};
		Boundary boundary = new BoundaryDeformablePeriodic(space, boxDim);
		Primitive primitive = new PrimitiveHexagonal(space, nC * aDim, nC * cDim);

		CoordinateDefinitionNitrogen coordinateDef = new CoordinateDefinitionNitrogen(this, box, primitive, basis, space);
		coordinateDef.setIsBetaHCP();
		coordinateDef.setOrientationVectorBeta(space);
		coordinateDef.initializeCoordinates(nCells);

		double[] u = new double[20];
		BetaPhaseLatticeParameter parameters = new BetaPhaseLatticeParameter();
		double[][] param = parameters.getParameter(density);

//			BetaPhaseLatticeParameterNA parameters = new BetaPhaseLatticeParameterNA();
//			double[][] param = parameters.getParameter(numMolecules);

		int kParam = 0;
		for (int i = 0; i < param.length; i++) {
			for (int j = 0; j < param[0].length; j++) {
				if (j < 3) {
					u[kParam] = 0.0;
				} else {
					u[kParam] = param[i][j];
				}
				kParam++;
			}
		}

		System.out.println("*************Parameters*************");
		for (int i = 0; i < u.length; i++) {
			System.out.print(u[i] + ", ");
			if ((i + 1) % 5 == 0) {
				System.out.println();
			}
		}

		int numDOF = coordinateDef.getCoordinateDim();
		double[] newU = new double[numDOF];
		if (true) {
			for (int j = 0; j < numDOF; j += 10) {
				if (j > 0 && j % (nC * 10) == 0) {
					j += nC * 10;
					if (j >= numDOF) {
						break;
					}
				}
				for (int k = 0; k < 10; k++) {
					newU[j + k] = u[k];
				}
			}

			for (int j = nC * 10; j < numDOF; j += 10) {
				if (j > nC * 10 && j % (nC * 10) == 0) {
					j += nC * 10;
					if (j >= numDOF) {
						break;
					}
				}
				for (int k = 0; k < 10; k++) {
					newU[j + k] = u[k + 10];
				}
			}
		}

		coordinateDef.setToU(box.getMoleculeList(), newU);
		coordinateDef.initNominalU(box.getMoleculeList());


		box.setBoundary(boundary);

		double rCScale = 0.475;
		double rc = aDim * nC * rCScale;
		System.out.println("Truncation Radius (" + rCScale + " Box Length): " + rc);
		P2Nitrogen potential = new P2Nitrogen(space, rc);
		potential.setBox(box);


		PotentialMasterListMolecular potentialMaster =
				new PotentialMasterListMolecular(this, rc, boxAgentSource, boxAgentManager, new NeighborListManagerSlantyMolecular.NeighborListSlantyAgentSourceMolecular(rc, space), space);
		potentialMaster.addPotential(potential, new ISpecies[]{species, species});

		MCMoveMoleculeCoupled move = new MCMoveMoleculeCoupled(potentialMaster, getRandom(), space);
		move.setBox(box);
		move.setPotential(potential);
		move.setDoExcludeNonNeighbors(true);

		IntegratorMC integrator = new IntegratorMC(potentialMaster, getRandom(), temperature, box);
		integrator.getMoveManager().addMCMove(move);

		int cellRange = 6;
		potentialMaster.setRange(rc);
		potentialMaster.setCellRange(cellRange);
		potentialMaster.getNeighborManager(box).reset();
		potential.setRange(Double.POSITIVE_INFINITY);

		int numNeigh = potentialMaster.getNeighborManager(box).getUpList(box.getMoleculeList().getMolecule(0))[0].getMoleculeCount();
		System.out.println("numNeigh: " + numNeigh);

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

        double temperature = 45; //in UNIT KELVIN
        double density = 0.023;
        double angle =1;
        long numSteps = 100000;
        int numMolecules = 1024;

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
	        
        System.out.println("Running beta-phase Nitrogen RP direct sampling simulation");
        System.out.println(numMolecules+" molecules at density "+density+" and temperature "+temperature + " K");
        System.out.println("perturbing from no rotational energy into rotational dof of angle= "+angle);
        System.out.println("with numStep of "+ numSteps);
        
        SimDirectBetaN2RPInitPert sim = new SimDirectBetaN2RPInitPert(Space.getInstance(3), numMolecules, density, Kelvin.UNIT.toSim(temperature), angle, numSteps);

        
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
