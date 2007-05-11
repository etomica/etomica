package etomica.normalmode;

import java.io.FileWriter;
import java.io.IOException;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomFactoryMono;
import etomica.atom.AtomType;
import etomica.atom.AtomTypeSphere;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.integrator.IntegratorHard;
import etomica.integrator.IntegratorMD;
import etomica.integrator.IntervalActionAdapter;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.lattice.crystal.PrimitiveFcc;
import etomica.phase.Phase;
import etomica.potential.P1HardPeriodic;
import etomica.potential.P2HardSphere;
import etomica.potential.Potential;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryDeformableLattice;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.IVector;
import etomica.space.Space;
import etomica.species.SpeciesSpheresMono;

/**
 * MD simulation of hard spheres in 1D or 3D with tabulation of the
 * collective-coordinate S-matrix. No graphic display of simulation.
 */
public class SimCalcS extends Simulation {

    public SimCalcS(Space space, int numAtoms, double density) {
        super(space, true, new PotentialMaster(space));

        defaults.makeLJDefaults();
        defaults.atomSize = 1.0;

        SpeciesSpheresMono species = new SpeciesSpheresMono(this);
        getSpeciesManager().addSpecies(species);

        phase = new Phase(this);
        addPhase(phase);
        phase.getAgent(species).setNMolecules(numAtoms);

        integrator = new IntegratorHard(potentialMaster, random, 0.04, 1.0);

        integrator.setIsothermal(false);
        activityIntegrate = new ActivityIntegrate(this, integrator);
        getController().addAction(activityIntegrate);
        // activityIntegrate.setMaxSteps(nSteps);

        Potential potential = new P2HardSphere(space, defaults.atomSize, false);
        AtomTypeSphere sphereType = (AtomTypeSphere) ((AtomFactoryMono) species
                .moleculeFactory()).getType();
        potentialMaster.addPotential(potential, new AtomType[] { sphereType,
                sphereType });

        int nCells;
        if (space.D() == 1) {
            primitive = new PrimitiveCubic(space, 1.0/density);
            nCells = numAtoms;
            bdry = new BoundaryRectangularPeriodic(space, getRandom(), numAtoms/density);
            ((IntegratorHard) integrator).setNullPotential(new P1HardPeriodic(space));
        } else {
            primitive = new PrimitiveFcc(space, 1);
            double v = primitive.unitCell().getVolume();
            primitive.scaleSize(Math.pow(v*density,-1.0/3.0));
            nCells = (int)Math.round(Math.pow(numAtoms, 1.0/3.0));
            bdry = new BoundaryDeformableLattice(primitive, getRandom(), new int[]{nCells,nCells,nCells});
        }
        phase.setBoundary(bdry);

        coordinateDefinition = new CoordinateDefinitionLeaf(phase, primitive);
        coordinateDefinition.initializeCoordinates(new int[]{nCells, nCells, nCells});
        
        integrator.setPhase(phase);
    }

    /**
     * @param args
     */
    public static void main(String[] args) {

        // defaults
        int D = 3;
        int nA = 27;
        double density = 1.3;
        if (D == 1) {
            nA = 3;
            density = 0.5;
        }

        double simTime = 100;

        // parse arguments
        if (args.length > 1) {
            density = Double.parseDouble(args[1]);
        }
        if (args.length > 2) {
            simTime = Double.parseDouble(args[2]);
        }
        if (args.length > 3) {
            nA = Integer.parseInt(args[3]);
        }
        String filename = "normal_modes" + D + "D_"+nA+"_"+((int)(density*100));
        if (args.length > 0) {
            filename = args[0];
        }

        System.out.println("Running "
                + (D == 1 ? "1D" : (D == 3 ? "FCC" : "2D hexagonal"))
                + " hard sphere simulation");
        System.out.println(nA + " atoms at density " + density);
        System.out.println(simTime + " time units");
        System.out.println("output data to " + filename);

        // construct simulation
        SimCalcS sim = new SimCalcS(Space.getInstance(D), nA, density);

        Primitive primitive = sim.primitive;

        // set up normal-mode meter
        MeterNormalMode meterNormalMode = new MeterNormalMode();
        meterNormalMode.setCoordinateDefinition(sim.coordinateDefinition);
        WaveVectorFactory waveVectorFactory;
        if (D == 1) {
            waveVectorFactory = new WaveVectorFactory1D();
        } else if (D == 2) {
            waveVectorFactory = null;
        } else {
            waveVectorFactory = new WaveVectorFactorySimple(primitive);
        }
        meterNormalMode.setWaveVectorFactory(waveVectorFactory);
        meterNormalMode.setPhase(sim.phase);

        IntervalActionAdapter fooAdapter = new IntervalActionAdapter(
                meterNormalMode);
        fooAdapter.setActionInterval(2);
        sim.integrator.addListener(fooAdapter);

        // MeterMomentumCOM meterCOM = new MeterMomentumCOM(sim.space);
        // MeterPositionCOM meterCOM = new MeterPositionCOM(sim.space);
        // DataSinkConsole console = new DataSinkConsole();
        // DataPump comPump = new DataPump(meterCOM,console);
        // IntervalActionAdapter comAdapter = new
        // IntervalActionAdapter(comPump);
        // sim.integrator.addListener(comAdapter);
        // meterCOM.setPhase(sim.phase);

        // start simulation
        int nSteps = (int) (simTime / sim.integrator.getTimeStep());
        sim.activityIntegrate.setMaxSteps(nSteps/10);
        sim.getController().actionPerformed();
        meterNormalMode.reset();
        sim.getController().reset();
        sim.activityIntegrate.setMaxSteps(nSteps);
        sim.getController().actionPerformed();

        // normalize averages
        DataGroup normalModeData = (DataGroup) meterNormalMode.getData();
        normalModeData.TE(1.0 / meterNormalMode.getCallCount());

        // write wave vectors (to filename.k) and simulation results (to
        // filename.S) to file
        IVector[] waveVectors = waveVectorFactory.getWaveVectors();
        double[] coefficients = waveVectorFactory.getCoefficients();

        try {
            int coordinateDim = meterNormalMode.getCoordinateDefinition()
                    .getCoordinateDim();
            FileWriter fileWriterK = new FileWriter(filename + ".k");
            FileWriter fileWriterS = new FileWriter(filename + ".S");
            for (int k = 0; k < waveVectors.length; k++) {
                // write the wavevector with its coefficient
                fileWriterK.write(Double.toString(coefficients[k]));
                for (int j = 0; j < waveVectors[k].getD(); j++) {
                    fileWriterK.write(" " + waveVectors[k].x(j));
                }
                fileWriterK.write("\n");
                if (D == 1) {
                    System.out.println(NormalModes1DHR.S1DHR(k + 1, nA / density, nA));
                }

                // write the (coordDim x coordDim) S array for the current
                // wavevector
                DataDoubleArray dataS = (DataDoubleArray) normalModeData.getData(k);
                for (int j = 0; j < coordinateDim; j++) {
                    fileWriterS.write(Double.toString(dataS.getValue(j * coordinateDim)));
                    for (int l = 1; l < coordinateDim; l++) {
                        fileWriterS.write(" " + dataS.getValue(j * coordinateDim + l));
                    }
                    fileWriterS.write("\n");
                }
            }
            fileWriterK.close();
            fileWriterS.close();
        } catch (IOException e) {
            throw new RuntimeException("Oops, failed to write data " + e);
        }

    }

    private static final long serialVersionUID = 1L;
    public IntegratorMD integrator;
    public ActivityIntegrate activityIntegrate;
    public Phase phase;
    public Boundary bdry;
    public Primitive primitive;
    public CoordinateDefinition coordinateDefinition;
}