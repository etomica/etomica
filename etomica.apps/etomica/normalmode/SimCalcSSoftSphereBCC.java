package etomica.normalmode;

import etomica.action.PDBWriter;
import etomica.action.activity.ActivityIntegrate;
import etomica.api.IAtomType;
import etomica.api.IBox;
import etomica.box.Box;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageCollapsing;
import etomica.data.DataPump;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.meter.MeterPressure;
import etomica.data.types.DataGroup;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveStepTracker;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisCubicBcc;
import etomica.lattice.crystal.BasisMonatomic;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.listener.IntegratorListenerAction;
import etomica.potential.P2SoftSphere;
import etomica.potential.P2SoftSphericalTruncatedShifted;
import etomica.potential.Potential2SoftSpherical;
import etomica.potential.PotentialMaster;
import etomica.potential.PotentialMasterMonatomic;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.species.SpeciesSpheresMono;

/**
 * MC simulation of BCC soft-sphere model in 3D with tabulation of the
 * collective-coordinate S-matrix. No graphic display of simulation.
 * 
 * @author Tai Boon Tan
 */
public class SimCalcSSoftSphereBCC extends Simulation {

    public SimCalcSSoftSphereBCC(Space _space, int numAtoms, double density, double temperature, int exponent) {
        super(_space, true);

        potentialMaster = new PotentialMasterMonatomic(this);

        SpeciesSpheresMono species = new SpeciesSpheresMono(this, space);
        getSpeciesManager().addSpecies(species);

        box = new Box(space);
        addBox(box);
        box.setNMolecules(species, numAtoms);

        integrator = new IntegratorMC(potentialMaster, getRandom(), temperature);
        MCMoveAtomCoupled move = new MCMoveAtomCoupled(potentialMaster, getRandom(), space);
        move.setStepSize(0.2);
        move.setStepSizeMax(0.5);
        integrator.getMoveManager().addMCMove(move);
        ((MCMoveStepTracker)move.getTracker()).setNoisyAdjustment(true);
        
        activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);
        // activityIntegrate.setMaxSteps(nSteps);

        if (space.D() == 1) {
            primitive = new PrimitiveCubic(space, 1.0/density);
            boundary = new BoundaryRectangularPeriodic(space, numAtoms/density);
            nCells = new int[]{numAtoms};
            basis = new BasisMonatomic(space);
        } else {
            double L = Math.pow(2.0/density, 1.0/3.0);
            primitive = new PrimitiveCubic(space, L);
            int n = (int)Math.round(Math.pow(numAtoms/2, 1.0/3.0));
            nCells = new int[]{n,n,n};
            boundary = new BoundaryRectangularPeriodic(space, n * L);
            basis = new BasisCubicBcc();
        }

        Potential2SoftSpherical potential = new P2SoftSphere(space, 1.0, 1.0, exponent);
        double truncationRadius = boundary.getBoxSize().getX(0) * 0.495;
        P2SoftSphericalTruncatedShifted pTruncated = new P2SoftSphericalTruncatedShifted(space, potential, truncationRadius);
        IAtomType sphereType = species.getLeafType();
        potentialMaster.addPotential(pTruncated, new IAtomType[] {sphereType, sphereType});
        move.setPotential(pTruncated);

        box.setBoundary(boundary);

        coordinateDefinition = new CoordinateDefinitionLeaf(this, box, primitive, basis, space);
        coordinateDefinition.initializeCoordinates(nCells);
        
        /*
         * 1-body Potential to Constraint the atom from moving too far 
         * 	away from its lattice-site
         */
       P1Constraint p1Constraint = new P1Constraint(space, primitive, box, coordinateDefinition);
       potentialMaster.addPotential(p1Constraint, new IAtomType[]{sphereType});
        
        integrator.setBox(box);
    }

    /**
     * @param args
     */
    public static void main(String[] args) {

        // defaults
        int D = 3;
        int nA = 128; 				//2*n^3 --- n be the number of unit cell
        double density = 2.206;		//set the density to melting density = 1.56*sqrt(2) Ref: Hoover(1971)
        double temperature = 0.1;	// varying temperature
        int exponent = 6;
        if (D == 1) {
            nA = 3;
            density = 0.5;
        }
        long simSteps = 1000000;

        // parse arguments
        if (args.length > 1) {
            density = Double.parseDouble(args[1]);
        }
        if (args.length > 2) {
            simSteps = Long.parseLong(args[2]);
        }
        if (args.length > 3) {
            nA = Integer.parseInt(args[3]);
        }
        if (args.length > 4) {
            temperature = Double.parseDouble(args[4]);
        }
        if (args.length > 5) {
        	exponent = Integer.parseInt(args[5]);
        }
        String filename = "CB_BCC_n"+exponent+"_T"+ (int)Math.round(temperature*10);
        if (args.length > 0) {
            filename = args[0];
        }

        System.out.println("Running "
                + (D == 1 ? "1D" : (D == 3 ? "BCC" : "2D hexagonal"))
                + " soft sphere simulation");
        System.out.println(nA + " atoms with exponent " + exponent+" and temperature "+temperature);
        System.out.println(simSteps+ " steps");
        System.out.println("output data to " + filename);

        // construct simulation
        SimCalcSSoftSphereBCC sim = new SimCalcSSoftSphereBCC(Space.getInstance(D), nA, density, temperature, exponent);

        // set up initial configuration and save nominal positions
        Primitive primitive = sim.primitive;
        
        /* 
         * set up normal-mode meter
         */
        
        
        MeterNormalMode meterNormalMode = new MeterNormalMode();
        meterNormalMode.setCoordinateDefinition(sim.coordinateDefinition);
        WaveVectorFactory waveVectorFactory;
        if (D == 1) {
            waveVectorFactory = new WaveVectorFactory1D();
        } else if (D == 2) {
            waveVectorFactory = null;
        } else {
            waveVectorFactory = new WaveVectorFactorySimple(primitive, sim.space);
        }
        
        meterNormalMode.setWaveVectorFactory(waveVectorFactory);
        meterNormalMode.setBox(sim.box);

        IntegratorListenerAction meterNormalModeListener = new IntegratorListenerAction(meterNormalMode);
        meterNormalModeListener.setInterval(nA);
        sim.integrator.getEventManager().addListener(meterNormalModeListener);
	
        MeterPressure meterPressure = new MeterPressure(sim.space);
        meterPressure.setIntegrator(sim.integrator);
        System.out.println("\nPressure Lattice: "+ meterPressure.getDataAsScalar());
        //System.exit(1);
        
        
        MeterPotentialEnergy meterEnergy = new MeterPotentialEnergy(sim.potentialMaster);
        meterEnergy.setBox(sim.box);
        System.out.println("Lattice Energy per particle: "+ meterEnergy.getDataAsScalar()/nA);
        System.out.println(" ");

        AccumulatorAverage pressureAverage = new AccumulatorAverageCollapsing();
	    DataPump pressurePump = new DataPump(meterPressure, pressureAverage);
	    
        IntegratorListenerAction pressurePumpListener = new IntegratorListenerAction(pressurePump);
        pressurePumpListener.setInterval(100);
	    sim.integrator.getEventManager().addListener(pressurePumpListener);
	    
        AccumulatorAverage energyAverage = new AccumulatorAverageCollapsing();
        DataPump energyPump = new DataPump(meterEnergy, energyAverage);
        
        IntegratorListenerAction energyPumpListener = new IntegratorListenerAction(energyPump);
        energyPumpListener.setInterval(100);
        sim.integrator.getEventManager().addListener(energyPumpListener);

      
        sim.activityIntegrate.setMaxSteps(simSteps/10);
        sim.getController().actionPerformed();
        System.out.println("equilibrated");
        
        long startTime = System.currentTimeMillis();
        System.out.println("\nStart Time: " +startTime);
        sim.integrator.getMoveManager().setEquilibrating(false);
        sim.getController().reset();
        meterNormalMode.reset();

        WriteS sWriter = new WriteS(sim.space);
        sWriter.setFilename(filename);
        sWriter.setOverwrite(true);
        sWriter.setMeter(meterNormalMode);
        sWriter.setWaveVectorFactory(waveVectorFactory);
        sWriter.setTemperature(temperature);
        IntegratorListenerAction sWriterListener = new IntegratorListenerAction(sWriter);
        sWriterListener.setInterval((int)simSteps/20);
        sim.integrator.getEventManager().addListener(sWriterListener);
        
        sim.activityIntegrate.setMaxSteps(simSteps);
        sim.getController().actionPerformed();
        PDBWriter pdbWriter = new PDBWriter(sim.box);
        pdbWriter.setFileName("calcS_n"+exponent+"_D"+density+".pdb");
        pdbWriter.actionPerformed();
       
        System.out.println("\nAverage Energy: "+ ((DataGroup)energyAverage.getData()).getValue(AccumulatorAverage.StatType.AVERAGE.index)+
        					" ,Error Energy: "+ ((DataGroup)energyAverage.getData()).getValue(AccumulatorAverage.StatType.ERROR.index));
        System.out.println(" ");
        
        System.out.println("Average-Pressure: "+ ((DataGroup)pressureAverage.getData()).getValue(AccumulatorAverage.StatType.AVERAGE.index) +
        					 " ,Error-Pressure: "+ ((DataGroup)pressureAverage.getData()).getValue(AccumulatorAverage.StatType.ERROR.index));
        
        long endTime = System.currentTimeMillis();
        System.out.println("End Time: " + endTime);
        System.out.println("Time taken: " + (endTime - startTime));
        
    }

    private static final long serialVersionUID = 1L;
    public IntegratorMC integrator;
    public ActivityIntegrate activityIntegrate;
    public IBox box;
    public Boundary boundary;
    public Primitive primitive;
    public Basis basis;
    public int[] nCells;
    public CoordinateDefinition coordinateDefinition;
    protected PotentialMaster potentialMaster;
}