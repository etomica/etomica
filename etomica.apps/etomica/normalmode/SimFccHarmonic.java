package etomica.normalmode;

import java.io.FileWriter;
import java.io.IOException;

import etomica.action.PhaseImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomFactoryMono;
import etomica.atom.AtomType;
import etomica.atom.AtomTypeSphere;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorAverage;
import etomica.data.DataFork;
import etomica.data.DataPump;
import etomica.data.AccumulatorAverage.StatType;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.integrator.IntegratorHard;
import etomica.integrator.IntegratorMD;
import etomica.integrator.IntervalActionAdapter;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.crystal.PrimitiveFcc;
import etomica.phase.Phase;
import etomica.potential.P2HardSphere;
import etomica.potential.Potential;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.species.SpeciesSpheresMono;

/**
 * Simulation to run sampling with the hard sphere potential, but measuring
 * the harmonic potential based on normal mode data from a previous simulation.
 * @author Andrew Schultz
 */
public class SimFccHarmonic extends Simulation {

    public SimFccHarmonic(Space space, int numAtoms, double density) {
        super(space, true, new PotentialMaster(space));

        defaults.makeLJDefaults();
        defaults.atomSize = 1.0;

        SpeciesSpheresMono species = new SpeciesSpheresMono(this);

        phase = new Phase(this);
        phase.getAgent(species).setNMolecules(numAtoms);

        integrator = new IntegratorHard(this);

        integrator.setIsothermal(false);
        activityIntegrate = new ActivityIntegrate(this,
                integrator);
        double timeStep = 0.04;
        integrator.setTimeStep(timeStep);
        getController().addAction(activityIntegrate);

        Potential potential = new P2HardSphere(space, defaults.atomSize, false);
        AtomTypeSphere sphereType = (AtomTypeSphere) ((AtomFactoryMono) species
                .moleculeFactory()).getType();
        potentialMaster.addPotential(potential, new AtomType[] { sphereType,
                sphereType });

        bdry = new BoundaryRectangularPeriodic(this);
        phase.setBoundary(bdry);
        phase.setDensity(density);

        PhaseImposePbc makeperiodic = new PhaseImposePbc(phase);
        integrator.addListener(makeperiodic);

        lattice = new LatticeCubicFcc();
        ConfigurationLattice config = new ConfigurationLattice(lattice);
        // config.setRescalingToFitVolume(false);

        config.initializeCoordinates(phase);

        // nan this section is a patch
        // first we find out the scaling used in
        // ConfigurationLattice/LatticeCubicFcc
        // then, we create a primitive fcc lattice, and scale it so we can use
        // it in pri.
        primitive = lattice.getPrimitiveFcc();
        ConfigurationLattice.MyLattice myLattice = (ConfigurationLattice.MyLattice) config
                .getLatticeMemento();
        Vector scaling = myLattice.latticeScaling;
        primitive.setCubicSize(primitive.getCubicSize()*scaling.x(0));

        // nan phase.setDensity(1.04);
        integrator.setPhase(phase);

    }

    /**
     * @param args
     */
    public static void main(String[] args) {
        int nA = 108;
        String filename = "normal_modes400";
        double density = 1.04;
        double simTime = 1000;
        if (args.length > 0) {
            filename = args[0];
        }
        if (args.length > 1) {
            density = Double.parseDouble(args[1]);
        }
        if (args.length > 2) {
            simTime = Double.parseDouble(args[2]);
        }
        if (args.length > 3) {
            nA = Integer.parseInt(args[3]);
        }
        
        System.out.println("Running FCC hard sphere simulation, measuring harmonic energy");
        System.out.println(nA+" atoms at density "+density);
        System.out.println(simTime+" time units");
        System.out.println("output data to "+filename);

        SimFccHarmonic sim = new SimFccHarmonic(Space3D.getInstance(), nA, density);
        
        double harmonicFudge = .25;
        
        double[][] omegaSquared = ArrayReader1D.getFromFile(filename+".val");
        for (int i=0; i<omegaSquared.length; i++) {
            for (int j=0; j<omegaSquared[i].length; j++) {
                // omega is sqrt(kT)/eigenvalue
                omegaSquared[i][j] = 1/omegaSquared[i][j]/harmonicFudge;
            }
        }
        double[][] waveVectorsAndCoefficients = ArrayReader1D.getFromFile(filename+".Q");
        Vector[] waveVectors = new Vector[waveVectorsAndCoefficients.length];
        double[] coefficients = new double[waveVectors.length];
        for (int i=0; i<waveVectors.length; i++) {
            coefficients[i] = waveVectorsAndCoefficients[i][0];
            waveVectors[i] = new Vector3D(waveVectorsAndCoefficients[i][1],
                    waveVectorsAndCoefficients[i][2],
                    waveVectorsAndCoefficients[i][3]);
        }
        double[][][] eigenvectors = ArrayReader2D.getFromFile(filename+".vec");

        MeterHarmonicEnergy harmonicEnergy = new MeterHarmonicEnergy();
        harmonicEnergy.setEigenvectors(eigenvectors);
        harmonicEnergy.setOmegaSquared(omegaSquared);
        harmonicEnergy.setWaveVectors(waveVectors, coefficients);
        harmonicEnergy.setNormalCoordWrapper(new NormalCoordLeaf(sim.getSpace()));
        harmonicEnergy.setPhase(sim.phase);
        DataFork harmonicFork = new DataFork();
        AccumulatorAverage harmonicAvg = new AccumulatorAverage(sim);
        DataPump pump = new DataPump(harmonicEnergy, harmonicFork);
        harmonicFork.addDataSink(harmonicAvg);
        IntervalActionAdapter adapter = new IntervalActionAdapter(pump);
        adapter.setActionInterval(2);
        sim.integrator.addListener(adapter);
        BoltzmannProcessor boltz = new BoltzmannProcessor();
        boltz.setTemperature(1.0);
        harmonicFork.addDataSink(boltz);
        AccumulatorAverage harmonicBoltzAvg = new AccumulatorAverage(50);
        boltz.setDataSink(harmonicBoltzAvg);
        
        MeterHarmonicSingleEnergy harmonicSingleEnergy = new MeterHarmonicSingleEnergy();
        harmonicSingleEnergy.setEigenvectors(eigenvectors);
        harmonicSingleEnergy.setOmegaSquared(omegaSquared);
        harmonicSingleEnergy.setWaveVectors(waveVectors, coefficients);
        harmonicSingleEnergy.setNormalCoordMapper(new NormalCoordLeaf(sim.getSpace()));
        harmonicSingleEnergy.setPhase(sim.phase);
        harmonicSingleEnergy.setTemperature(1.0);
        AccumulatorAverage harmonicSingleAvg = new AccumulatorAverage(sim);
        pump = new DataPump(harmonicSingleEnergy, harmonicSingleAvg);
        adapter = new IntervalActionAdapter(pump);
        adapter.setActionInterval(2);
        sim.integrator.addListener(adapter);

        int nSteps = (int) (simTime / sim.integrator.getTimeStep());

        sim.activityIntegrate.setMaxSteps(nSteps);
        
        sim.getController().actionPerformed();

        double avgHarmonicEnergy = ((DataDouble)((DataGroup)harmonicAvg.getData()).getData(AccumulatorAverage.StatType.AVERAGE.index)).x;
        double errorHarmonicEnergy = ((DataDouble)((DataGroup)harmonicAvg.getData()).getData(AccumulatorAverage.StatType.ERROR.index)).x;
        System.out.println("avg harmonic energy: "+avgHarmonicEnergy+" +/- "+errorHarmonicEnergy);
        
        DataDoubleArray harmonicModesAvg = (DataDoubleArray)((DataGroup)harmonicSingleAvg.getData()).getData(StatType.AVERAGE.index);
        DataDoubleArray harmonicModesErr = (DataDoubleArray)((DataGroup)harmonicSingleAvg.getData()).getData(StatType.ERROR.index);

        double deltaA = 0;
        double deltaAerr = 0;
        int nData = harmonicModesAvg.getLength();
        for (int i=0; i<nData; i++) {
            deltaA += Math.log(harmonicModesAvg.getValue(i));
            deltaAerr += harmonicModesErr.getValue(i)/harmonicModesAvg.getValue(i);
        }
        
        System.out.println("Harmonic free energy correction (independent approx): "+deltaA+" +/- "+deltaAerr);

        deltaA = ((DataDouble)((DataGroup)harmonicBoltzAvg.getData()).getData(StatType.AVERAGE.index)).x;
        deltaAerr = ((DataDouble)((DataGroup)harmonicBoltzAvg.getData()).getData(StatType.ERROR.index)).x/deltaA;
        deltaA = Math.log(deltaA);
        
        System.out.println("Harmonic free energy correction: "+deltaA+" +/- "+deltaAerr);

        try {
            // energy -- e? u?
            FileWriter fileWriterE = new FileWriter(filename+".e");
            int[] idx = new int[2];
            for (int i=0; i<harmonicModesAvg.getArrayShape(0); i++) {
                idx[0] = i;
                idx[1] = 0;
                fileWriterE.write(Double.toString(harmonicModesAvg.getValue(idx)));
                for (int j=1; j<harmonicModesAvg.getArrayShape(1); j++) {
                    idx[1] = j;
                    fileWriterE.write(" "+harmonicModesAvg.getValue(idx));
                }
                fileWriterE.write("\n");
            }
            fileWriterE.close();
        }
        catch (IOException e) {
            throw new RuntimeException("Oops, failed to write data "+e);
        }
    }

    private static final long serialVersionUID = 1L;
    public IntegratorMD integrator;
    public ActivityIntegrate activityIntegrate;
    public Phase phase;
    public BoundaryRectangularPeriodic bdry;
    public LatticeCubicFcc lattice;
    public PrimitiveFcc primitive;
}
