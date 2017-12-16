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
 * @author Tai Boon Tan
 */
public class SimDirectBetaN2RPAngleToNoAngle extends Simulation {

    public SimDirectBetaN2RPAngleToNoAngle(Space space, int numMolecules, double density, double temperature, double angle, long numSteps) {
        super(space);
        
        BoxAgentSourceCellManagerListMolecular boxAgentSourceTarg = new BoxAgentSourceCellManagerListMolecular(this, null, space);
        BoxAgentManager<NeighborCellManagerMolecular> boxAgentManagerTarg = new BoxAgentManager<NeighborCellManagerMolecular>(boxAgentSourceTarg, NeighborCellManagerMolecular.class, this);
 
        BoxAgentSourceCellManagerListMolecular boxAgentSourceRef = new BoxAgentSourceCellManagerListMolecular(this, null, space);
        BoxAgentManager<NeighborCellManagerMolecular> boxAgentManagerRef = new BoxAgentManager<NeighborCellManagerMolecular>(boxAgentSourceRef, NeighborCellManagerMolecular.class, this);
        
        SpeciesN2 species = new SpeciesN2(space);
		addSpecies(species);

        // TARGET
        boxTarg = new Box(space);
        addBox(boxTarg);
        boxTarg.setNMolecules(species, numMolecules);

    	double ratio = 1.631;
		double aDim = Math.pow(4.0/(Math.sqrt(3.0)*ratio*density), 1.0/3.0);
		double cDim = aDim*ratio;
		System.out.println("\naDim: " + aDim + " ;cDim: " + cDim);
		int nC = (int)Math.pow(numMolecules/1.999999999, 1.0/3.0);
		
		Basis basisHCP = new BasisHcp();
		BasisBigCell basis = new BasisBigCell(space, basisHCP, new int[]{nC,nC,nC});
        
		Vector[] boxDim = new Vector[3];
		boxDim[0] = space.makeVector(new double[]{nC*aDim, 0, 0});
		boxDim[1] = space.makeVector(new double[]{-nC*aDim*Math.cos(Degree.UNIT.toSim(60)), nC*aDim*Math.sin(Degree.UNIT.toSim(60)), 0});
		boxDim[2] = space.makeVector(new double[]{0, 0, nC*cDim});
		
		int[] nCells = new int[]{1,1,1};
		Boundary boundary = new BoundaryDeformablePeriodic(space, boxDim);
		Primitive primitive = new PrimitiveHexagonal(space, nC*aDim, nC*cDim);
		
		CoordinateDefinitionNitrogen coordinateDefTarg = new CoordinateDefinitionNitrogen(this, boxTarg, primitive, basis, space);
		coordinateDefTarg.setIsBetaHCP();
		coordinateDefTarg.setOrientationVectorBeta(space);
		coordinateDefTarg.initializeCoordinates(nCells);
		
	    boxTarg.setBoundary(boundary);
	    
		double rCScale = 0.475;
		double rc = aDim*nC*rCScale;
		System.out.println("Truncation Radius (" + rCScale +" Box Length): " + rc);
		P2Nitrogen potentialTarg = new P2Nitrogen(space, rc);
		potentialTarg.setBox(boxTarg);

		PRotConstraint pRotConstraintTarg = new PRotConstraint(space,coordinateDefTarg, boxTarg);
		pRotConstraintTarg.setConstraintAngle(angle);
		
		PotentialMasterListMolecular potentialMasterTarg = 
			new PotentialMasterListMolecular(this, rc, boxAgentSourceTarg, boxAgentManagerTarg, new NeighborListManagerSlantyMolecular.NeighborListSlantyAgentSourceMolecular(rc, space), space);
		potentialMasterTarg.addPotential(potentialTarg, new ISpecies[]{species, species});
		//potentialMasterTarg.addPotential(pRotConstraintTarg,new ISpecies[]{species} );
	       
		MCMoveMoleculeCoupledInitPert moveTarg = new MCMoveMoleculeCoupledInitPert(potentialMasterTarg, getRandom(),space, coordinateDefTarg);
		moveTarg.setBox(boxTarg);
		moveTarg.setPotential(potentialTarg);
		moveTarg.setDoExcludeNonNeighbors(true);
		
//		MCMoveRotateMolecule3D rotateTarg = new MCMoveRotateMolecule3D(potentialMasterTarg, getRandom(), space);
//		rotateTarg.setBox(boxTarg);
		MCMoveRotateMolecule3DFixedAngle rotateFixedAngle
		= new MCMoveRotateMolecule3DFixedAngle(potentialMasterTarg, getRandom(), space, angle, coordinateDefTarg, boxTarg);
		rotateFixedAngle.setBox(boxTarg);
		
        IntegratorMC integratorTarg = new IntegratorMC(potentialMasterTarg, getRandom(), temperature);
        integratorTarg.getMoveManager().addMCMove(moveTarg);
		//integratorTarg.getMoveManager().addMCMove(rotateTarg);
        integratorTarg.getMoveManager().addMCMove(rotateFixedAngle);
	    
		integratorTarg.setBox(boxTarg);

        // Reference System
		PotentialMasterListMolecular potentialMasterRef = 
			new PotentialMasterListMolecular(this, rc, boxAgentSourceRef, boxAgentManagerRef, new NeighborListManagerSlantyMolecular.NeighborListSlantyAgentSourceMolecular(rc, space), space);
		potentialMasterRef.addPotential(potentialTarg, new ISpecies[]{species, species});
		potentialMasterRef.addPotential(pRotConstraintTarg,new ISpecies[]{species} );
		
		int cellRange = 6;
	    potentialMasterTarg.setRange(rc);
	    potentialMasterTarg.setCellRange(cellRange); 
	    potentialMasterTarg.getNeighborManager(boxTarg).reset();

	        
//	    int potentialCells = potentialMasterTarg.getNbrCellManager(boxTarg).getLattice().getSize()[0];
//	    if (potentialCells < cellRange*2+1) {
//	      throw new RuntimeException("oops ("+potentialCells+" < "+(cellRange*2+1)+")");
//	    }
	    potentialMasterRef.setRange(rc);
	    potentialMasterRef.setCellRange(cellRange); 
	    potentialMasterRef.getNeighborManager(boxTarg).reset();
	    potentialTarg.setRange(Double.POSITIVE_INFINITY); 
	    
	    MeterPotentialEnergy meterPERef = new MeterPotentialEnergy(potentialMasterRef);
        meterPERef.setBox(boxTarg);
        double latticeEnergy = meterPERef.getDataAsScalar();
        System.out.println("lattice energy (sim unit): " + latticeEnergy);
		System.out.println("lattice energy per mol(sim unit): " + latticeEnergy/numMolecules);
		
		MeterBoltzmannDirect meterBoltzmann = new MeterBoltzmannDirect(integratorTarg, meterPERef);
		
//        MeterRotPerturbMolecule meterBoltzmannRotPerb = new MeterRotPerturbMolecule(integratorTarg, potentialMasterRef, species, space, this, coordinateDefTarg);
//        meterBoltzmannRotPerb.setLatticeEnergy(latticeEnergy);
		int numBlock = 100;
		int interval = numMolecules;
		long blockSize = numSteps/(numBlock*interval);
        if (blockSize == 0) blockSize = 1;
        System.out.println("block size "+blockSize+" interval "+interval);
        boltzmannAverage = new AccumulatorAverageFixed(blockSize);
        
        DataPump boltzmannPump = new DataPump(meterBoltzmann, boltzmannAverage);
        IntegratorListenerAction boltzmannPumpListener = new IntegratorListenerAction(boltzmannPump, 100);
        integratorTarg.getEventManager().addListener(boltzmannPumpListener);

        activityIntegrate = new ActivityIntegrate(integratorTarg);
        getController().addAction(activityIntegrate);
    }

    public void initializeConfigFromFile(String fname){
        ConfigurationFile config = new ConfigurationFile(fname);
        config.initializeCoordinates(boxTarg);
    }
    
    public void writeConfiguration(String fname){
        WriteConfiguration writeConfig = new WriteConfiguration(space);
        writeConfig.setBox(boxTarg);
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
        double density = 0.025;
        double angle =1.0;
        long numSteps = 1000000;
        int numMolecules = 432;

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
		String configFileName = "configT"+temperature+"_Angle"+angle;
	        
        System.out.println("Running beta-phase Nitrogen RP direct sampling simulation");
        System.out.println(numMolecules+" molecules at density "+density+" and temperature "+temperature + " K");
        System.out.println("perturbing from no rotational energy into rotational dof of angle= "+angle);
        System.out.println("with numStep of "+ numSteps);
        
        SimDirectBetaN2RPAngleToNoAngle sim = new SimDirectBetaN2RPAngleToNoAngle(Space.getInstance(3), numMolecules, density, Kelvin.UNIT.toSim(temperature), angle, numSteps);

        
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
    protected Box boxTarg;

}
