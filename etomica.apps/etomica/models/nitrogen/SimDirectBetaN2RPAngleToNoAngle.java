package etomica.models.nitrogen;

import etomica.action.activity.ActivityIntegrate;
import etomica.api.ISpecies;
import etomica.api.IVector;
import etomica.box.Box;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.DataPump;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.types.DataGroup;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveRotateMolecule3D;
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
import etomica.util.ParameterBase;
import etomica.util.ReadParameters;

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
        Box boxTarg = new Box(space);
        addBox(boxTarg);
        boxTarg.setNMolecules(species, numMolecules);

    	double ratio = 1.631;
		double aDim = Math.pow(4.0/(Math.sqrt(3.0)*ratio*density), 1.0/3.0);
		double cDim = aDim*ratio;
		System.out.println("\n\naDim: " + aDim + " ;cDim: " + cDim);
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
		coordinateDefTarg.setIsBeta();
		coordinateDefTarg.setOrientationVectorBeta(space);
		coordinateDefTarg.initializeCoordinates(nCells);
		
	    boxTarg.setBoundary(boundary);
	    
		double rCScale = 0.485;
		double rc = aDim*nC*rCScale;
		System.out.println("Truncation Radius (" + rCScale +" Box Length): " + rc);
		P2Nitrogen potentialTarg = new P2Nitrogen(space, rc);
		potentialTarg.setBox(boxTarg);

		PRotConstraint pRotConstraintTarg = new PRotConstraint(space,coordinateDefTarg, boxTarg);
		pRotConstraintTarg.setConstraintAngle(angle);
		
		potentialMasterTarg.addPotential(potentialTarg, new ISpecies[]{species, species});
		potentialMasterTarg.addPotential(pRotConstraintTarg,new ISpecies[]{species} );
		
		MCMoveMoleculeCoupled moveTarg = new MCMoveMoleculeCoupled(potentialMasterTarg, getRandom(),space);
		moveTarg.setBox(boxTarg);
		moveTarg.setPotential(potentialTarg);
		
		MCMoveRotateMolecule3D rotateTarg = new MCMoveRotateMolecule3D(potentialMasterTarg, getRandom(), space);
		rotateTarg.setBox(boxTarg);
		
        IntegratorMC integratorTarg = new IntegratorMC(potentialMasterTarg, getRandom(), temperature);
        integratorTarg.getMoveManager().addMCMove(moveTarg);
		integratorTarg.getMoveManager().addMCMove(rotateTarg);
	    
		integratorTarg.setBox(boxTarg);
		
	    MeterPotentialEnergy meterPETarg = new MeterPotentialEnergy(potentialMasterTarg);
        meterPETarg.setBox(boxTarg);
        double latticeEnergy = meterPETarg.getDataAsScalar();
        System.out.println("lattice energy (sim unit): " + latticeEnergy);
        
        // Reference System
		potentialMasterRef.addPotential(potentialTarg, new ISpecies[]{species, species});
        MeterRotPerturbMolecule meterBoltzmannRotPerb = new MeterRotPerturbMolecule(integratorTarg, potentialMasterRef, species, space, this, coordinateDefTarg);
        meterBoltzmannRotPerb.setLatticeEnergy(latticeEnergy);
        boltzmannAverage = new AccumulatorAverageFixed(100);
        
        DataPump boltzmannPump = new DataPump(meterBoltzmannRotPerb, boltzmannAverage);
        IntegratorListenerAction boltzmannPumpListener = new IntegratorListenerAction(boltzmannPump, 100);
        integratorTarg.getEventManager().addListener(boltzmannPumpListener);

        activityIntegrate = new ActivityIntegrate(integratorTarg);
        getController().addAction(activityIntegrate);
    }

    /**
     * @param args filename containing simulation parameters
     * @see SimDirectBetaN2RPAngleToNoAngle.SimOverlapParam
     */
    public static void main(String[] args) {
        //set up simulation parameters
        SimOverlapParam params = new SimOverlapParam();
        
        String inputFilename = null;
        if (args.length > 0) {
            inputFilename = args[0];
        }
        if (inputFilename != null) {
            ReadParameters readParameters = new ReadParameters(inputFilename, params);
            readParameters.readParameters();
        }
        
        double density = params.density;
        double angle = params.angle;
        long numSteps = params.numSteps;
        int numMolecules = params.numMolecules;
        double temperature = params.temperature;
  
        System.out.println("Running beta-phase Nitrogen RP direct sampling simulation");
        System.out.println(numMolecules+" molecules at density "+density+" and temperature "+temperature + " K");
        System.out.print("perturbing from angle=" + angle + " into no rotational d.o.f.");
        
        SimDirectBetaN2RPAngleToNoAngle sim = new SimDirectBetaN2RPAngleToNoAngle(Space.getInstance(3), numMolecules, density, Kelvin.UNIT.toSim(temperature), angle);

        sim.activityIntegrate.setMaxSteps(numSteps/10);
        sim.getController().actionPerformed();     
        System.out.println("equilibration finished");
        sim.getController().reset();
     
 
        long startTime = System.currentTimeMillis();
        System.out.println("Start Time: " + startTime);
       
        sim.activityIntegrate.setMaxSteps(numSteps);
        sim.getController().actionPerformed();

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
    
    /**
     * Inner class for parameters understood by the SimOverlapBetaN2RP constructor
     */
    public static class SimOverlapParam extends ParameterBase {
        public int numMolecules = 432;
        public double density = 0.025;
        public double angle = 1;
        public int D = 3;
        public long numSteps =50000;
        public double temperature = 40;
    }
}
