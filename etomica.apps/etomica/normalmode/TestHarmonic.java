package etomica.normalmode;

import java.io.FileWriter;
import java.io.IOException;

import etomica.action.PhaseImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorAverage;
import etomica.data.DataFork;
import etomica.data.DataHistogram;
import etomica.data.DataPump;
import etomica.data.AccumulatorAverage.StatType;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.graphics.DisplayBoxesCAE;
import etomica.graphics.DisplayPlot;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.integrator.IntervalActionAdapter;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.crystal.PrimitiveFcc;
import etomica.phase.Phase;
import etomica.potential.P2HardSphere;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.Species;
import etomica.species.SpeciesSpheresMono;
import etomica.util.DoubleRange;
import etomica.util.HistogramSimple;

/**
 * Simulation to sample harmonic potential
 */
public class TestHarmonic extends Simulation {

    public TestHarmonic(Space space, int numAtoms, String filename) {
        super(space, true, new PotentialMaster(space));

        defaults.makeLJDefaults();
        defaults.atomSize = 1.0;

        species = new SpeciesSpheresMono(this);

        phase = new Phase(this);
        phase.getAgent(species).setNMolecules(numAtoms);

        integrator = new IntegratorMC(this);

        activityIntegrate = new ActivityIntegrate(this,
                integrator);
        getController().addAction(activityIntegrate);

        MCMoveHarmonic move = new MCMoveHarmonic(potentialMaster);
        integrator.getMoveManager().addMCMove(move);
        
        bdry = new BoundaryRectangularPeriodic(this);
        phase.setBoundary(bdry);
        phase.setDensity(1.04);

        PhaseImposePbc makeperiodic = new PhaseImposePbc(phase);
        integrator.addListener(makeperiodic);

        double[][] eigenValues = ArrayReader1D.getFromFile(filename+".val");
        Vector[] q = ArrayReader1D.getVectorsFromFile(filename+".Q");
        double[][][] eigenvectors = ArrayReader2D.getFromFile(filename+".vec");
        
        double harmonicFudge = 1;
        for (int i=0; i<eigenValues.length; i++) {
            for (int j=0; j<eigenValues[i].length; j++) {
                eigenValues[i][j] = harmonicFudge*eigenValues[i][j];
            }
        }
        
        move.setEigenValues(eigenValues);
        move.setEigenVectors(eigenvectors);
        move.setWaveVectors(q);
        
        lattice = new LatticeCubicFcc();
        config = new ConfigurationLattice(lattice);
        // config.setRescalingToFitVolume(false);

        config.initializeCoordinates(phase);

        move.setPhase(phase);
        
        integrator.setPhase(phase);
    }

    /**
     * @param args
     */
    public static void main(String[] args) {
        int nA = 108;
        boolean graphic = true;
        String filename = "normal_modes400";
        if (args.length > 0) {
            filename = args[0];
        }
        TestHarmonic sim = new TestHarmonic(Space3D.getInstance(), nA, filename);
        
        P2HardSphere p2HardSphere = new P2HardSphere(sim.space, 1.0, true);
        sim.potentialMaster.addPotential(p2HardSphere, new AtomType[]{sim.species.getMoleculeType(),sim.species.getMoleculeType()});
        
        MeterPotentialEnergy meterPE = new MeterPotentialEnergy(sim.potentialMaster);
        meterPE.setPhase(sim.phase);
        BoltzmannProcessor bp = new BoltzmannProcessor();
        bp.setTemperature(1);
        DataPump pump = new DataPump(meterPE,bp);
        AccumulatorAverage avgBoltzmann = new AccumulatorAverage(1);
        bp.setDataSink(avgBoltzmann);
        avgBoltzmann.setPushInterval(5);
        IntervalActionAdapter iaa = new IntervalActionAdapter(pump);
        sim.integrator.addListener(iaa);
        
        NormalCoordLeaf normalCoordLeaf = new NormalCoordLeaf(sim.space);
        double harmonicFudge = 1;
        Vector[] q = ArrayReader1D.getVectorsFromFile(filename+".Q");
        double[][][] eigenvectors = ArrayReader2D.getFromFile(filename+".vec");
        double[][] omegaSquared = ArrayReader1D.getFromFile(filename+".val");
        for (int i=0; i<omegaSquared.length; i++) {
            for (int j=0; j<omegaSquared[i].length; j++) {
                // omega is sqrt(kT)/eigenvalue
                omegaSquared[i][j] = 1/omegaSquared[i][j]/harmonicFudge;
            }
        }
        MeterHarmonicEnergy harmonicEnergy = new MeterHarmonicEnergy();
        harmonicEnergy.setEigenvectors(eigenvectors);
        harmonicEnergy.setOmegaSquared(omegaSquared);
        harmonicEnergy.setWaveVectors(q);
        harmonicEnergy.setPhase(sim.phase);
        harmonicEnergy.setNormalCoordWrapper(normalCoordLeaf);
        DataFork harmonicFork = new DataFork();
        AccumulatorAverage harmonicAvg = new AccumulatorAverage(5);
        pump = new DataPump(harmonicEnergy, harmonicFork);
        harmonicFork.addDataSink(harmonicAvg);
        IntervalActionAdapter adapter = new IntervalActionAdapter(pump);
        adapter.setActionInterval(2);
        sim.integrator.addListener(adapter);

        MeterHarmonicSingleEnergy harmonicSingleEnergy = new MeterHarmonicSingleEnergy();
        harmonicSingleEnergy.setEigenvectors(eigenvectors);
        harmonicSingleEnergy.setOmegaSquared(omegaSquared);
        harmonicSingleEnergy.setWaveVectors(q);
        harmonicSingleEnergy.setPhase(sim.phase);
        harmonicSingleEnergy.setTemperature(1.0);
        harmonicSingleEnergy.setNormalCoordMapper(normalCoordLeaf);
//        DataProcessorFunction harmonicLog = new DataProcessorFunction(new Function.Log());
        AccumulatorAverage harmonicSingleAvg = new AccumulatorAverage(5);
        DataHistogram harmonicSingleHistogram = new DataHistogram(new HistogramSimple.Factory(50, new DoubleRange(0, 1)));
        pump = new DataPump(harmonicSingleEnergy, harmonicSingleHistogram);
//        harmonicLog.setDataSink(harmonicSingleHistogram);
        harmonicSingleHistogram.setDataSink(harmonicSingleAvg);
        iaa= new IntervalActionAdapter(pump);
        iaa.setActionInterval(1);
        sim.integrator.addListener(iaa);

        if(graphic){
            SimulationGraphic simG = new SimulationGraphic(sim);
            DisplayBoxesCAE boxesPE = new DisplayBoxesCAE();
            boxesPE.setAccumulator(avgBoltzmann);
            boxesPE.setPrecision(6);
            simG.add(boxesPE);

            DisplayBoxesCAE harmonicBoxes = new DisplayBoxesCAE();
            harmonicBoxes.setAccumulator(harmonicAvg);
            simG.add(harmonicBoxes);

            DisplayPlot harmonicPlot = new DisplayPlot();
            harmonicPlot.setDoLegend(false);
            harmonicSingleAvg.addDataSink(harmonicPlot.getDataSet().makeDataSink(), new StatType[]{StatType.AVERAGE});
            simG.add(harmonicPlot);

            simG.makeAndDisplayFrame();
        } else {
            PrimitiveFcc primitive = sim.lattice.getPrimitiveFcc();
            ConfigurationLattice.MyLattice myLattice = (ConfigurationLattice.MyLattice) sim.config
                    .getLatticeMemento();
            Vector scaling = myLattice.latticeScaling;
            primitive.setCubicSize(primitive.getCubicSize()*scaling.x(0));

            MeterNormalMode meterNormalMode = new MeterNormalMode();
            meterNormalMode.setPhase(sim.phase);
            meterNormalMode.setNormalCoordWrapper(normalCoordLeaf);
            meterNormalMode.setWaveVectorFactory(new WaveVectorFactoryFcc(primitive));

            int nSteps = 100;

            sim.activityIntegrate.setMaxSteps(nSteps);

            String outFilename = "normal_modes_harmonic";
            if (args.length > 0) {
                filename = args[0];
            }
            
            sim.activityIntegrate.setMaxSteps(nSteps);

            if (true) {
                IntervalActionAdapter fooAdapter = new IntervalActionAdapter(meterNormalMode);
                fooAdapter.setActionInterval(1);
                sim.integrator.addListener(fooAdapter);
            }
            
            sim.getController().actionPerformed();
            
            if (true) {
                DataGroup normalModeData = (DataGroup)meterNormalMode.getData();
                normalModeData.TE(1.0/(sim.phase.getSpeciesMaster().moleculeCount()*meterNormalMode.getCallCount()));
                int normalDim = meterNormalMode.getNormalCoordWrapper().getNormalDim();
                
                Vector[] waveVectors = meterNormalMode.getWaveVectors();
                
                try {
                    FileWriter fileWriterQ = new FileWriter(outFilename+".Q");
                    FileWriter fileWriterS = new FileWriter(outFilename+".S");
                    for (int i=0; i<waveVectors.length; i++) {
                        fileWriterQ.write(Double.toString(waveVectors[i].x(0)));
                        for (int j=1; j<waveVectors[i].D(); j++) {
                            fileWriterQ.write(" "+waveVectors[i].x(j));
                        }
                        fileWriterQ.write("\n");
                        DataDoubleArray dataS = (DataDoubleArray)normalModeData.getData(i);
                        for (int k=0; k<normalDim; k++) {
                            fileWriterS.write(Double.toString(dataS.getValue(k*normalDim)));
                            for (int l=1; l<normalDim; l++) {
                                fileWriterS.write(" "+dataS.getValue(k*normalDim+l));
                            }
                            fileWriterS.write("\n");
                        }
                    }
                    fileWriterQ.close();
                    fileWriterS.close();
                }
                catch (IOException e) {
                    throw new RuntimeException("Oops, failed to write data "+e);
                }
            }
        }

    }

    private static final long serialVersionUID = 1L;
    public IntegratorMC integrator;
    public ActivityIntegrate activityIntegrate;
    public Phase phase;
    public BoundaryRectangularPeriodic bdry;
    public LatticeCubicFcc lattice;
    public ConfigurationLattice config;
    public Species species;
}