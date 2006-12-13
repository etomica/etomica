package etomica.normalmode;

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
import etomica.data.types.DataDouble;
import etomica.data.types.DataGroup;
import etomica.graphics.DisplayBoxesCAE;
import etomica.graphics.DisplayPlot;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.integrator.IntervalActionAdapter;
import etomica.lattice.BravaisLattice;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.LatticeCubicSimple;
import etomica.phase.Phase;
import etomica.potential.P2HardSphere;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.Species;
import etomica.species.SpeciesSpheresMono;
import etomica.util.DoubleRange;
import etomica.util.HistogramSimple;

/**
 * Simulation to sample harmonic potential
 */
public class TestHarmonic extends Simulation {

    public TestHarmonic(Space space, int numAtoms, double density, String filename) {
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
        phase.setDensity(density);

        PhaseImposePbc makeperiodic = new PhaseImposePbc(phase);
        integrator.addListener(makeperiodic);

        double harmonicFudge = .25;
        double[][] eigenValues = ArrayReader1D.getFromFile(filename+".val");
        for (int i=0; i<eigenValues.length; i++) {
            for (int j=0; j<eigenValues[i].length; j++) {
                // omega is sqrt(kT)/eigenvalue
                eigenValues[i][j] = eigenValues[i][j]*harmonicFudge;
            }
        }
        double[][] waveVectorsAndCoefficients = ArrayReader1D.getFromFile(filename+".Q");
        Vector[] waveVectors = new Vector[waveVectorsAndCoefficients.length];
        double[] coefficients = new double[waveVectors.length];
        double[] justWaveVector = new double[space.D()];
        for (int i=0; i<waveVectors.length; i++) {
            coefficients[i] = waveVectorsAndCoefficients[i][0];
            for (int j=0; j<space.D(); j++) {
                justWaveVector[j] = waveVectorsAndCoefficients[i][j+1];
            }
            waveVectors[i] = Space.makeVector(justWaveVector); 
        }
        double[][][] eigenvectors = ArrayReader2D.getFromFile(filename+".vec");
        
        move.setEigenValues(eigenValues);
        move.setEigenVectors(eigenvectors);
        move.setWaveVectors(waveVectors);
        move.setWaveVectorCoefficients(coefficients);
        move.setNormalCoordWrapper(new NormalCoordLeaf(space));
        
        if (space.D() == 1) {
            lattice = new LatticeCubicSimple(1,phase.getBoundary().getDimensions().x(0)/numAtoms);
        }
        else {
            lattice = new LatticeCubicFcc();
        }
        config = new ConfigurationLattice(lattice);

        config.initializeCoordinates(phase);

        move.setPhase(phase);
        
        integrator.setPhase(phase);
    }

    /**
     * @param args
     */
    public static void main(String[] args) {
        int D = 1;
        int nA = 108;
        double density = 1.04;
        if (D == 1) {
            nA = 5;
            density = 0.5;
        }
        boolean graphic = false;
        String filename = "normal_modes400";
        if (args.length > 0) {
            filename = args[0];
        }
        TestHarmonic sim = new TestHarmonic(Space.getInstance(D), nA, density, filename);
        
        P2HardSphere p2HardSphere = new P2HardSphere(sim.getSpace(), 1.0, true);
        sim.getPotentialMaster().addPotential(p2HardSphere, new AtomType[]{sim.species.getMoleculeType(),sim.species.getMoleculeType()});
        
        MeterPotentialEnergy meterPE = new MeterPotentialEnergy(sim.getPotentialMaster());
        meterPE.setPhase(sim.phase);
        BoltzmannProcessor bp = new BoltzmannProcessor();
        bp.setTemperature(1);
        DataPump pump = new DataPump(meterPE,bp);
        AccumulatorAverage avgBoltzmann = new AccumulatorAverage(1);
        bp.setDataSink(avgBoltzmann);
        avgBoltzmann.setPushInterval(5);
        IntervalActionAdapter iaa = new IntervalActionAdapter(pump);
        sim.integrator.addListener(iaa);
        
        NormalCoordLeaf normalCoordLeaf = new NormalCoordLeaf(sim.getSpace());
        double harmonicFudge = .25;
        double[][] waveVectorsAndCoefficients = ArrayReader1D.getFromFile(filename+".Q");
        Vector[] waveVectors = new Vector[waveVectorsAndCoefficients.length];
        double[] coefficients = new double[waveVectors.length];
        double[] justWaveVector = new double[D];
        for (int i=0; i<waveVectors.length; i++) {
            coefficients[i] = waveVectorsAndCoefficients[i][0];
            for (int j=0; j<D; j++) {
                justWaveVector[j] = waveVectorsAndCoefficients[i][j+1];
            }
            waveVectors[i] = Space.makeVector(justWaveVector); 
        }
        double[][][] eigenvectors = ArrayReader2D.getFromFile(filename+".vec");
        //these are actually eigenvalues
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
        harmonicEnergy.setWaveVectors(waveVectors, coefficients);
        harmonicEnergy.setNormalCoordWrapper(normalCoordLeaf);
        harmonicEnergy.setPhase(sim.phase);
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
        harmonicSingleEnergy.setWaveVectors(waveVectors, coefficients);
        harmonicSingleEnergy.setNormalCoordMapper(normalCoordLeaf);
        harmonicSingleEnergy.setTemperature(1.0);
        harmonicSingleEnergy.setPhase(sim.phase);
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

            int nSteps = 1000000;

            sim.activityIntegrate.setMaxSteps(nSteps);
            
            sim.getController().actionPerformed();
            
            DataGroup boltzmannData = (DataGroup)avgBoltzmann.getData();
            double pNotOverlap = ((DataDouble)boltzmannData.getData(StatType.AVERAGE.index)).x;
            double pError = ((DataDouble)boltzmannData.getData(StatType.ERROR.index)).x;
            
            System.out.println("avg HS Boltzmann factor "+pNotOverlap+" +/- "+pError);
            
            System.out.println("free energy contribution "+(-Math.log(pNotOverlap))+" +/- "+(pError/pNotOverlap));
        }

    }

    private static final long serialVersionUID = 1L;
    public IntegratorMC integrator;
    public ActivityIntegrate activityIntegrate;
    public Phase phase;
    public BoundaryRectangularPeriodic bdry;
    public BravaisLattice lattice;
    public ConfigurationLattice config;
    public Species species;
}