package etomica.models.nitrogen;

import java.io.File;

import etomica.action.WriteConfiguration;
import etomica.action.activity.ActivityIntegrate;
import etomica.api.IBox;
import etomica.api.ISpecies;
import etomica.api.IVector;
import etomica.box.Box;
import etomica.config.ConfigurationFile;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.DataPump;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.types.DataGroup;
import etomica.integrator.IntegratorMC;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisHcp;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveHexagonal;
import etomica.listener.IntegratorListenerAction;
import etomica.normalmode.BasisBigCell;
import etomica.normalmode.MCMoveMoleculeCoupled;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryDeformablePeriodic;
import etomica.space.Space;
import etomica.units.Degree;
import etomica.units.Kelvin;

/**
 * Direct Sampling for Rotational Perturbation
 *  from constraint angle to no rotational d.o.f.
 *   
 * @author Tai Boon Tan
 */
public class SimDirectBetaN2RPAngleToNoAngle extends Simulation {

    public SimDirectBetaN2RPAngleToNoAngle(Space space, int numMolecules, double density, double temperature, double angle) {
        super(space);
        
        PotentialMaster potentialMasterTarg = new PotentialMaster();
        PotentialMaster potentialMasterRef = new PotentialMaster();
        
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
        
		IVector[] boxDim = new IVector[3];
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
		
		potentialMasterTarg.addPotential(potentialTarg, new ISpecies[]{species, species});
		//potentialMasterTarg.addPotential(pRotConstraintTarg,new ISpecies[]{species} );
		
		MCMoveMoleculeCoupled moveTarg = new MCMoveMoleculeCoupled(potentialMasterTarg, getRandom(),space);
		moveTarg.setBox(boxTarg);
		moveTarg.setPotential(potentialTarg);
		
//		MCMoveRotateMolecule3D rotateTarg = new MCMoveRotateMolecule3D(potentialMasterTarg, getRandom(), space);
//		rotateTarg.setBox(boxTarg);
		
		MCMoveRotateMolecule3DConstraint rotateConst
		= new MCMoveRotateMolecule3DConstraint(potentialMasterTarg, getRandom(), space, angle, coordinateDefTarg, boxTarg);
		rotateConst.setBox(boxTarg);
		
        IntegratorMC integratorTarg = new IntegratorMC(potentialMasterTarg, getRandom(), temperature);
        integratorTarg.getMoveManager().addMCMove(moveTarg);
		//integratorTarg.getMoveManager().addMCMove(rotateTarg);
        integratorTarg.getMoveManager().addMCMove(rotateConst);
	    
		integratorTarg.setBox(boxTarg);

        
        // Reference System
		potentialMasterRef.addPotential(potentialTarg, new ISpecies[]{species, species});
		potentialMasterRef.addPotential(pRotConstraintTarg,new ISpecies[]{species} );
		
	    MeterPotentialEnergy meterPERef = new MeterPotentialEnergy(potentialMasterRef);
        meterPERef.setBox(boxTarg);
        double latticeEnergy = meterPERef.getDataAsScalar();
        System.out.println("lattice energy per mol(Kelvin): " + Kelvin.UNIT.fromSim(latticeEnergy)/numMolecules);
        System.out.println("lattice energy per mol(sim unit): " + latticeEnergy/numMolecules);
		
		MeterBoltzmannDirect meterBoltzmann = new MeterBoltzmannDirect(integratorTarg, meterPERef);
		
//        MeterRotPerturbMolecule meterBoltzmannRotPerb = new MeterRotPerturbMolecule(integratorTarg, potentialMasterRef, species, space, this, coordinateDefTarg);
//        meterBoltzmannRotPerb.setLatticeEnergy(latticeEnergy);
        boltzmannAverage = new AccumulatorAverageFixed(100);
        
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
        double angle = 1.0;
        long numSteps = 100000;
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
        System.out.println("perturbing from angle =" + angle + " into no rotational d.o.f.");
        System.out.println("with numStep of "+ numSteps);
        
        SimDirectBetaN2RPAngleToNoAngle sim = new SimDirectBetaN2RPAngleToNoAngle(Space.getInstance(3), numMolecules, density, Kelvin.UNIT.toSim(temperature), angle);

        
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
        double average = ((DataGroup)sim.boltzmannAverage.getData()).getValue(AccumulatorAverage.StatType.AVERAGE.index);
        double error = ((DataGroup)sim.boltzmannAverage.getData()).getValue(AccumulatorAverage.StatType.ERROR.index);
        
        System.out.println("boltzmann average: " + average + " ;err: " + error);
        
        long endTime = System.currentTimeMillis();
        System.out.println("\nEnd Time: " + endTime);
        System.out.println("Time taken (s): " + (endTime - startTime)/1000);
       
    }

    private static final long serialVersionUID = 1L;
    protected ActivityIntegrate activityIntegrate;
    protected AccumulatorAverageFixed boltzmannAverage;
    protected IBox boxTarg;

}
