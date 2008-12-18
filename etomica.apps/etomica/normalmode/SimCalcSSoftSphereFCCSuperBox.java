package etomica.normalmode;

import etomica.action.PDBWriter;
import etomica.action.activity.ActivityIntegrate;
import etomica.api.IAtomTypeLeaf;
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
import etomica.lattice.crystal.BasisCubicFcc;
import etomica.lattice.crystal.BasisMonatomic;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.potential.P2SoftSphere;
import etomica.potential.P2SoftSphericalTruncatedShifted;
import etomica.potential.Potential2SoftSpherical;
import etomica.potential.PotentialMasterMonatomic;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.species.SpeciesSpheresMono;

/**
 * MC simulation of FCC soft-sphere model in 3D with tabulation of the
 * collective-coordinate S-matrix. No graphic display of simulation.
 * 
 * with SuperBox
 * 
 * @author Tai Boon Tan
 */
public class SimCalcSSoftSphereFCCSuperBox extends Simulation {

    public SimCalcSSoftSphereFCCSuperBox(Space _space, int numAtoms, double density, double temperature, int exponent) {
        super(_space, true);


        potentialMaster = new PotentialMasterMonatomic(this);

        SpeciesSpheresMono species = new SpeciesSpheresMono(this, space);
        getSpeciesManager().addSpecies(species);

        box = new Box(space);
        addBox(box);
        box.setNMolecules(species, numAtoms);

        if (space.D() == 1) {
            primitive = new PrimitiveCubic(space, 1.0/density);
            boundary = new BoundaryRectangularPeriodic(space, numAtoms/density);
            nCells = new int[]{numAtoms};
            basis = new BasisMonatomic(space);
        } else {
            double L = Math.pow(4.0/density, 1.0/3.0);
            primitive = new PrimitiveCubic(space, L);
            int n = (int)Math.round(Math.pow(numAtoms/4, 1.0/3.0));
            nCells = new int[]{n,n,n};
            boundary = new BoundaryRectangularPeriodic(space, n * L);
            basis = new BasisCubicFcc();
        }

        Potential2SoftSpherical potential = new P2SoftSphere(space, 1.0, 1.0, exponent);
        
        double truncationRadius = boundary.getDimensions().x(0) * 0.495;
        P2SoftSphericalTruncatedShifted pTruncated = new P2SoftSphericalTruncatedShifted(space, potential, truncationRadius);
        //potentialMaster.lrcMaster().setEnabled(false); //turn off the long-range correction ::updated 7/4/2008 
        
        IAtomTypeLeaf sphereType = species.getLeafType();
        potentialMaster.addPotential(pTruncated, new IAtomTypeLeaf[] {sphereType, sphereType});
        
        box.setBoundary(boundary);

        coordinateDefinition = new CoordinateDefinitionLeaf(this, box, primitive, basis, space);
        coordinateDefinition.initializeCoordinates(nCells);
        
        integrator = new IntegratorMC(potentialMaster, getRandom(), temperature);
        MCMoveAtomSuperBox move = new MCMoveAtomSuperBox(potentialMaster, getRandom(), space, coordinateDefinition);
        move.setStepSize(0.2);
        move.setStepSizeMax(0.5);
        integrator.getMoveManager().addMCMove(move);
        ((MCMoveStepTracker)move.getTracker()).setNoisyAdjustment(true);
        
        activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);
        
        move.setPotential(pTruncated);
        // activityIntegrate.setMaxSteps(nSteps);
        
        /*
         * 1-body Potential to Constraint the atom from moving too far 
         * 	away from its lattice-site
         */
       P1Constraint p1Constraint = new P1Constraint(space, primitive, box, coordinateDefinition);
       potentialMaster.addPotential(p1Constraint, new IAtomTypeLeaf[]{sphereType});
      
       integrator.setBox(box);
    }

    /**
     * @param args
     */
    public static void main(String[] args) {

        // defaults
        int D = 3;
        int nA = 32;
        double density = 1.256;
        double temperature = 1.2;
        int exponent = 12;
        if (D == 1) {
            nA = 3;
            density = 1.0;
        }
        long simSteps =1000000;

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
        String filename = "FCC_SoftSphere_n" + exponent + "_T"+temperature;
        if (args.length > 0) {
            filename = args[0];
        }

        System.out.println("Running "
                + (D == 1 ? "1D" : (D == 3 ? "FCC" : "2D hexagonal"))
                + " soft sphere simulation");
        System.out.println(nA + " atoms with exponent " + exponent+" and density "+density);
        System.out.println("isotherm temperature at "+temperature);
        System.out.println(simSteps+ " steps");
        System.out.println("output data to " + filename);

        // construct simulation
        SimCalcSSoftSphereFCCSuperBox sim = new SimCalcSSoftSphereFCCSuperBox(Space.getInstance(D), nA*27, density, temperature, exponent);

        // set up initial configuration and save nominal positions
        Primitive primitive = sim.primitive;

        // set up normal-mode meter
        MeterNormalMode meterNormalMode = new MeterNormalMode();
        meterNormalMode.setCoordinateDefinition(sim.coordinateDefinition);
        WaveVectorFactory waveVectorFactory;
        if (D == 1) {
            waveVectorFactory = new WaveVectorFactory1D(D);
        } else if (D == 2) {
            waveVectorFactory = null;
        } else {
            waveVectorFactory = new WaveVectorFactorySuperBox(primitive, sim.space);
        }
        meterNormalMode.setWaveVectorFactory(waveVectorFactory);
        meterNormalMode.setBox(sim.box);

        sim.integrator.addIntervalAction(meterNormalMode);
        sim.integrator.setActionInterval(meterNormalMode, nA);

        MeterPressure meterPressure = new MeterPressure(sim.space);
        meterPressure.setIntegrator(sim.integrator);
        System.out.println("\nPressure Lattice: "+ meterPressure.getDataAsScalar());
        //System.exit(1);
        
        
        MeterPotentialEnergy meterEnergy = new MeterPotentialEnergy(sim.potentialMaster);
        meterEnergy.setBox(sim.box);
        System.out.println("Lattice Energy per particle: "+ meterEnergy.getDataAsScalar()/(nA*27));
        System.out.println(" ");
        //System.exit(1);
        
        AccumulatorAverage pressureAverage = new AccumulatorAverageCollapsing();
	    DataPump pressurePump = new DataPump(meterPressure, pressureAverage);
	    
	    sim.integrator.addIntervalAction(pressurePump);
	    sim.integrator.setActionInterval(pressurePump, 100);
	    
        AccumulatorAverage energyAverage = new AccumulatorAverageCollapsing();
        DataPump energyPump = new DataPump(meterEnergy, energyAverage);
        
        sim.integrator.addIntervalAction(energyPump);
        sim.integrator.setActionInterval(energyPump, 100);
        
        sim.activityIntegrate.setMaxSteps(simSteps/10);  //simSteps/10
        sim.getController().actionPerformed();
        System.out.println("equilibrated");
        sim.integrator.getMoveManager().setEquilibrating(false);
        sim.getController().reset();
        meterNormalMode.reset();

        WriteS sWriter = new WriteS(sim.space);
        sWriter.setFilename(filename);
        sWriter.setOverwrite(true);
        sWriter.setMeter(meterNormalMode);
        sWriter.setWaveVectorFactory(waveVectorFactory);
        sWriter.setTemperature(temperature);
        sim.integrator.addIntervalAction(sWriter);
        sim.integrator.setActionInterval(sWriter, (int)simSteps/10);
        
        sim.activityIntegrate.setMaxSteps(simSteps);
        sim.getController().actionPerformed();
        PDBWriter pdbWriter = new PDBWriter(sim.box);
        pdbWriter.setFileName("calcS_n"+exponent+"_T"+temperature+".pdb");
        pdbWriter.actionPerformed();
        
        System.out.println("Average Energy: "+ ((DataGroup)energyAverage.getData()).getValue(AccumulatorAverage.StatType.AVERAGE.index)/27);
        System.out.println("Error Energy: "+ ((DataGroup)energyAverage.getData()).getValue(AccumulatorAverage.StatType.ERROR.index)/27);
        System.out.println(" ");
        
        System.out.println("Average Pressure: "+ ((DataGroup)pressureAverage.getData()).getValue(AccumulatorAverage.StatType.AVERAGE.index));
        System.out.println("Error Pressure: "+ ((DataGroup)pressureAverage.getData()).getValue(AccumulatorAverage.StatType.ERROR.index));
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
    public PotentialMasterMonatomic potentialMaster;
}