package etomica.starpolymer;


import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.DiameterHashByType;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.history.HistoryCollapsingAverage;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.meter.MeterRadiusGyration;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayPlot;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorMD;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.*;
import etomica.potential.compute.PotentialCompute;
import etomica.potential.compute.PotentialComputeAggregate;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularNonperiodic;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.species.SpeciesGeneral;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.util.random.RandomMersenneTwister;

import java.util.ArrayList;
import java.util.List;

public class StarPolymerMD extends Simulation {

    public Box box;
    public int f, l;
    public SpeciesGeneral species;
    public IntegratorVelocityVerlet integratorMD;
    public P2Fene potentialFene;
    public P2WCA potentialWCA;
    public PotentialCompute potentialMaster;

    public StarPolymerMD(int f, int l, double temperature, double tStep, boolean useNbrs) {
        super(Space3D.getInstance());
        setRandom(new RandomMersenneTwister(2));
        this.f = f;
        this.l = l;
        species = SpeciesPolymerMono.create(getSpace(), AtomType.simpleFromSim(this), f, l)
                .setDynamic(true)
                .build();
        addSpecies(species);

        box = this.makeBox(new BoundaryRectangularNonperiodic(space));
        box.getBoundary().setBoxSize(new Vector3D(1.5 * l * 2, 1.5 * l * 2, 1.5 * l * 2));
        box.setNMolecules(species, 1);
        List<int[]> bondedPairs = getPairArray(f, l);
        PotentialMasterBonding pmBonding = new PotentialMasterBonding(getSpeciesManager(), box);
        PotentialMaster potentialMasterPair;
        if (useNbrs) {
            potentialMasterPair = new PotentialMasterList(getSpeciesManager(), box, 1, 2, pmBonding.getBondingInfo());
        } else {
            potentialMasterPair = new PotentialMaster(getSpeciesManager(), box, pmBonding.getBondingInfo());
        }
        potentialMaster = new PotentialComputeAggregate(pmBonding, potentialMasterPair);
        integratorMD = new IntegratorVelocityVerlet(potentialMaster, this.getRandom(), tStep, temperature, box);
        integratorMD.setIsothermal(true);
        integratorMD.setThermostat(IntegratorMD.ThermostatType.ANDERSEN);
        integratorMD.setThermostatInterval(500);

        this.getController().addActivity(new ActivityIntegrate(integratorMD));

        potentialFene = new P2Fene();
        potentialWCA = new P2WCA();

        P2SoftSphericalSum pBonding = new P2SoftSphericalSum(potentialFene, potentialWCA);

        pmBonding.setBondingPotentialPair(species, pBonding, bondedPairs);

        AtomType type = species.getAtomType(0);
        potentialMasterPair.setPairPotential(type, type, potentialWCA);

    }

    public static void main(String[] args) {
        StarPolymerMD.PolymerParam params = new StarPolymerMD.PolymerParam();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        } else {
            params.fp = "./resource/rg.dat";
            params.numSteps = 1000000L;
            params.fv = 3;
            params.lv = 34;
            params.xsteps = 1000000L;
        }

        String fp = params.fp;
        long steps = params.numSteps;
        int f = params.fv;
        int l = params.lv;
        long xsteps = params.xsteps;

        double temperature = 1.0;
        final int dataInterval = 10;
        boolean graphics = true;
        double tStep = 0.005;
        boolean useNbrs = true;

        final StarPolymerMD sim = new StarPolymerMD(f, l, temperature, tStep, useNbrs);

        final MeterPotentialEnergyFromIntegrator meterPE = new MeterPotentialEnergyFromIntegrator(sim.integratorMD);
        double u = meterPE.getDataAsScalar();
        System.out.println("Initial Potential energy: " + u);

        DataFork forkPE = new DataFork();
        AccumulatorAverageCollapsing accPE = new AccumulatorAverageCollapsing();
        AccumulatorHistory histPE = new AccumulatorHistory(new HistoryCollapsingAverage());
        forkPE.addDataSink(accPE);
        forkPE.addDataSink(histPE);
        DataPumpListener pumpPF = new DataPumpListener(meterPE, forkPE, dataInterval);
        sim.integratorMD.getEventManager().addListener(pumpPF);
        sim.integratorMD.reset();

        MeterRadiusGyration meterRG = new MeterRadiusGyration(sim.getBox(0));
        System.out.println(String.format("Squared Radius of Gyration: %.8f\n", meterRG.getDataAsScalar()));

        DataFork rgFork = new DataFork();
        AccumulatorAverageCollapsing accRG = new AccumulatorAverageCollapsing();
        AccumulatorHistory histRG = new AccumulatorHistory(new HistoryCollapsingAverage());
        rgFork.addDataSink(accRG);
        rgFork.addDataSink(histRG);
        DataPumpListener rgPump = new DataPumpListener(meterRG, rgFork, dataInterval);
        sim.integratorMD.getEventManager().addListener(rgPump);
        sim.integratorMD.reset();

        MeterBondLength meterBL = new MeterBondLength(sim.box, f, l);

        if (graphics) {
            final String APP_NAME = "StarPolymer-MD";
//            sim.getBox(0).getBoundary().setBoxSize(Vector.of(new double[]{50, 50, 50}));
            final SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, APP_NAME, 1);
            DisplayBox displayBox = simGraphic.getDisplayBox(sim.getBox(0));
            displayBox.setShowBoundary(false);
            simGraphic.getController().getReinitButton().setPostAction(simGraphic.getPaintAction(sim.box));
            simGraphic.makeAndDisplayFrame(APP_NAME);
            DiameterHashByType dh = (DiameterHashByType) simGraphic.getDisplayBox(sim.box).getDiameterHash();
            dh.setDiameter(sim.species.getAtomType(0), 1.0);
//            sim.ai.setMaxSteps(1);

            DisplayPlot pePlot = new DisplayPlot();
            histPE.addDataSink(pePlot.getDataSet().makeDataSink());
            pePlot.setLabel("PE");
            simGraphic.add(pePlot);

            sim.integratorMD.getEventManager().addListener(new IntegratorListenerAction(meterBL, dataInterval));
            DisplayPlot blPlot = new DisplayPlot();
            blPlot.getPlot().setYLog(true);
            DataPumpListener blPump = new DataPumpListener(meterBL, blPlot.getDataSet().makeDataSink(), dataInterval);
            sim.integratorMD.getEventManager().addListener(blPump);
            blPlot.setLabel("bond length");
            simGraphic.add(blPlot);

            DisplayPlot rgPlot = new DisplayPlot();
            histRG.addDataSink(rgPlot.getDataSet().makeDataSink());
            rgPlot.setLabel("Radius of Gyration^2");
            simGraphic.add(rgPlot);

            return;
        }

        System.out.println("MD starts...");
        long t1 = System.currentTimeMillis();
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integratorMD, steps));
        long t2 = System.currentTimeMillis();
        System.out.println("MD finished! ");
        System.out.println("time : " + (t2 - t1) / 1000.0);

        if (true) {
            System.out.println();
            System.out.println("MD begin for dumping rg^2 ");
            int interval = 1000;
            DataLogger dataLogger = new DataLogger();
            DataPumpListener rgPumpListener = new DataPumpListener(meterRG, dataLogger, interval);
            sim.integratorMD.getEventManager().addListener(rgPumpListener);
            dataLogger.setFileName(fp);
            dataLogger.setAppending(false);
            DataArrayWriter writer = new DataArrayWriter();
            writer.setIncludeHeader(false);
            dataLogger.setDataSink(writer);
            long t3 = System.currentTimeMillis();
            sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integratorMD, xsteps));
            dataLogger.cleanUp();
            System.out.println("Dumping finished! ");
            System.out.println("time: " + (t3 - t1) / 1000.0);
        }
        System.out.println();
        System.out.println("Final Potential energy: " + meterPE.getDataAsScalar());
        System.out.println(String.format("Final rg^2 : %.8f\n", meterRG.getDataAsScalar()));

    }

    private ArrayList<int[]> getPairArray(int f, int l) {
        ArrayList<int[]> bondedPairArray = new ArrayList<>();

        for (int k = 0; k < f; k++) {
            int[] temp = new int[]{0, 1+k*l};
            bondedPairArray.add(temp);
        }
        for (int k = 0; k < f; k++) {
            for (int i=0; i<l-1; i++) {
                int a1 = 1 + k*l + i;
                bondedPairArray.add(new int[]{a1, a1+1});
            }
        }
        return bondedPairArray;
    }

    public static class PolymerParam extends ParameterBase {
        public long numSteps = 1000000;
        public String fp = null;
        public int fv = 5;
        public int lv = 40;
        public int size = 1000;
        public int interval = 1000;
        public long xsteps = 1000000;
    }

}
