package etomica.normalmode;

import java.io.FileWriter;
import java.io.IOException;

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
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.graphics.DisplayBoxesCAE;
import etomica.graphics.DisplayPlot;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.integrator.IntervalActionAdapter;
import etomica.lattice.BravaisLattice;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.LatticeCubicSimple;
import etomica.lattice.crystal.Primitive;
import etomica.phase.Phase;
import etomica.potential.P2HardSphere;
import etomica.potential.P2XOrder;
import etomica.potential.Potential;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.IVector;
import etomica.space.Space;
import etomica.species.Species;
import etomica.species.SpeciesSpheresMono;
import etomica.util.DoubleRange;
import etomica.util.HistogramSimple;

/**
 * Simulation to sample harmonic potential
 */
public class SimHarmonic extends Simulation {

    public SimHarmonic(Space space, int numAtoms, double density, String filename, double harmonicFudge) {
        super(space, true, new PotentialMaster(space));

        int D = space.D();
        
        defaults.makeLJDefaults();
        defaults.atomSize = 1.0;

        species = new SpeciesSpheresMono(this);
        getSpeciesManager().addSpecies(species);

        phase = new Phase(this);
        phase.getAgent(species).setNMolecules(numAtoms);

        integrator = new IntegratorMC(this);

        activityIntegrate = new ActivityIntegrate(this, integrator);
        getController().addAction(activityIntegrate);

        MCMoveHarmonic move = new MCMoveHarmonic(potentialMaster, getRandom());
        integrator.getMoveManager().addMCMove(move);
        
        bdry = new BoundaryRectangularPeriodic(this);
        phase.setBoundary(bdry);
        phase.setDensity(density);

        if (space.D() == 1) {
            lattice = new LatticeCubicSimple(1,phase.getBoundary().getDimensions().x(0)/numAtoms);
        }
        else {
            lattice = new LatticeCubicFcc();
        }
        config = new ConfigurationLattice(lattice);

        config.initializeCoordinates(phase);
        
        if(D == 1) {
            normalModes = new NormalModes1DHR();
        } else {
            normalModes = new NormalModesFromFile(filename, D);
        }
        
        WaveVectorFactory waveVectorFactory = normalModes.getWaveVectorFactory();
        waveVectorFactory.makeWaveVectors(phase);
        move.setEigenValues(normalModes.getEigenvalues(phase), waveVectorFactory.getCoefficients());
        move.setEigenVectors(normalModes.getEigenvectors(phase));
        move.setWaveVectors(waveVectorFactory.getWaveVectors());
        move.setWaveVectorCoefficients(waveVectorFactory.getCoefficients());
        move.setCoordinateDefinition(new CoordinateDefinitionLeaf(space));
        
        move.setPhase(phase);
        
        integrator.setPhase(phase);
    }

    /**
     * @param args
     */
    public static void main(String[] args) {
        
        //set up simulation parameters
        int D = 1;
        int nA = 108;
        double density = 1.04;
        if (D == 1) {
            nA = 3;
            density = 0.5;
        }
        boolean graphic = true;
        String filename = "normal_modes1D";
        if (args.length > 0) {
            filename = args[0];
        }
        double harmonicFudge = 1;
        
        //construct simulation
        SimHarmonic sim = new SimHarmonic(Space.getInstance(D), nA, density, filename, harmonicFudge);
        
        //add hard potentials for FEP calculations.  With de novo sampling potential is not otherwise used.
        P2HardSphere p2HardSphere = new P2HardSphere(sim.getSpace(), 1.0, true);
        sim.getPotentialMaster().addPotential(p2HardSphere, new AtomType[]{sim.species.getMoleculeType(),sim.species.getMoleculeType()});
        if (sim.getSpace().D() == 1) {
            Potential pOrder = new P2XOrder(sim.getSpace());
            sim.potentialMaster.addPotential(pOrder, new AtomType[] {sim.species.getMoleculeType(),sim.species.getMoleculeType()});
        }
        
        //meters for FEP calculations
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

//         MeterMomentumCOM meterCOM = new MeterMomentumCOM(sim.space);
//         MeterPositionCOM meterCOM = new MeterPositionCOM(sim.space);
//         DataSinkConsole console = new DataSinkConsole();
//         DataProcessorFunction filter = new DataProcessorFunction(new Function.Chop());
//         DataPump comPump = new DataPump(meterCOM,filter);
//         filter.setDataSink(console);
//         IntervalActionAdapter comAdapter = new IntervalActionAdapter(comPump);
//         sim.integrator.addListener(comAdapter);
//         meterCOM.setPhase(sim.phase);

        //set up things for determining energy of harmonic system
        //read and set up wave vectors
        CoordinateDefinitionLeaf coordinateDefinitionLeaf = new CoordinateDefinitionLeaf(sim.getSpace());

        //meter for harmonic system energy, sent to direct and to boltzmann average
        MeterHarmonicEnergy harmonicEnergy = new MeterHarmonicEnergy(coordinateDefinitionLeaf, sim.normalModes);
        harmonicEnergy.setPhase(sim.phase);
        DataFork harmonicFork = new DataFork();
        AccumulatorAverage harmonicAvg = new AccumulatorAverage(5);
        pump = new DataPump(harmonicEnergy, harmonicFork);
        harmonicFork.addDataSink(harmonicAvg);
        IntervalActionAdapter adapter = new IntervalActionAdapter(pump);
        adapter.setActionInterval(1);
        sim.integrator.addListener(adapter);

        //histogram energy of individual modes
        MeterHarmonicSingleEnergy harmonicSingleEnergy = new MeterHarmonicSingleEnergy(coordinateDefinitionLeaf, sim.normalModes);
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
        
        //set up measurement of S matrix, to check that configurations are generated as expected
        Primitive primitive = sim.lattice.getPrimitive();
        if (D == 3) {
            primitive = ((LatticeCubicFcc)sim.lattice).getPrimitiveFcc();
        }
        ConfigurationLattice.MyLattice myLattice = (ConfigurationLattice.MyLattice) sim.config.getLatticeMemento();
        IVector scaling = myLattice.latticeScaling;
        primitive.scaleSize(scaling.x(0));
        MeterNormalMode meterNormalMode = new MeterNormalMode();
        meterNormalMode.setCoordinateDefinition(coordinateDefinitionLeaf);
        WaveVectorFactory waveVectorFactory = sim.normalModes.getWaveVectorFactory();
        meterNormalMode.setWaveVectorFactory(waveVectorFactory);
        meterNormalMode.setPhase(sim.phase);


        if(graphic){
            //graphic simulation -- set up window
//            sim.getDefaults().pixelUnit = new Pixel(0.05);
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
            //not graphic, so run simulation batch
            //S data is written to file
            int nSteps = 100000;
            sim.activityIntegrate.setMaxSteps(nSteps);

            IntervalActionAdapter fooAdapter = new IntervalActionAdapter(meterNormalMode);
            fooAdapter.setActionInterval(2);
            sim.integrator.addListener(fooAdapter);
            sim.getController().actionPerformed();
            
            DataGroup normalModeData = (DataGroup)meterNormalMode.getData();
            normalModeData.TE(1.0/(sim.phase.getSpeciesMaster().moleculeCount()*meterNormalMode.getCallCount()));
            int normalDim = meterNormalMode.getCoordinateDefinition().getCoordinateDim();
            
            
            try {
                FileWriter fileWriterS = new FileWriter(filename+"_redo.S");
                for (int i=0; i<sim.normalModes.getWaveVectorFactory().getWaveVectors().length; i++) {
                    DataDoubleArray dataS = (DataDoubleArray)normalModeData.getData(i);
                    for (int k=0; k<normalDim; k++) {
                        fileWriterS.write(Double.toString(dataS.getValue(k*normalDim)));
                        for (int l=1; l<normalDim; l++) {
                            fileWriterS.write(" "+dataS.getValue(k*normalDim+l));
                        }
                        fileWriterS.write("\n");
                    }
                }
                fileWriterS.close();
            }
            catch (IOException e) {
                throw new RuntimeException("Oops, failed to write data "+e);
            }

            
            DataGroup boltzmannData = (DataGroup)avgBoltzmann.getData();
            double pNotOverlap = ((DataDouble)boltzmannData.getData(StatType.AVERAGE.index)).x;
            double pError = ((DataDouble)boltzmannData.getData(StatType.ERROR.index)).x;
            
            System.out.println("avg HS Boltzmann factor "+pNotOverlap+" +/- "+pError);
            
            System.out.println("free energy contribution "+(-Math.log(pNotOverlap))+" +/- "+(pError/pNotOverlap));
            System.out.println("free energy contribution per molecule "+(-Math.log(pNotOverlap)/nA)+" +/- "+(pError/pNotOverlap)/nA);
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
    public NormalModes normalModes;
}