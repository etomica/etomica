package etomica.starpolymer;

import etomica.action.IAction;
import etomica.action.WriteConfiguration;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.DiameterHashByType;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.history.HistoryCollapsingAverage;
import etomica.data.meter.MeterPotentialEnergyFromIntegratorFasterer;
import etomica.data.meter.MeterRadiusGyration;
import etomica.data.types.DataGroup;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayPlot;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorMCFasterer;
import etomica.nbr.cell.PotentialMasterCellFasterer;
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

import java.io.File;
import java.util.ArrayList;
import java.util.List;

public class StarPolymerMCFasterer extends Simulation {

    public Box box;
    public int f, l;
    public SpeciesGeneral species;
    public IntegratorMCFasterer integratorMC;
    public P2Fene potentialFene;
    public P2WCA potentialWCA;
    public PotentialCompute potentialMaster;

    public StarPolymerMCFasterer(int f, int l, double temperature) {
        super(Space3D.getInstance());
        this.f = f;
        this.l = l;
        species = SpeciesPolymerMono.create(this.getSpace(), AtomType.simpleFromSim(this), f, l)
                .setDynamic(true)
                .build();
        addSpecies(species);

        box = this.makeBox(new BoundaryRectangularNonperiodic(space));
        box.getBoundary().setBoxSize(new Vector3D(1.5 * l * 2, 1.5 * l * 2, 1.5 * l * 2));
        box.setNMolecules(species, 1);
        List<int[]> bondedPairs = getPairArray(f, l);

        boolean useNbrs = false;
        if (l < 30) useNbrs = false;
        PotentialMasterBonding pmBonding = new PotentialMasterBonding(getSpeciesManager(), box);
        PotentialMasterFasterer potentialMasterPair;
        if (useNbrs) {
            potentialMasterPair = new PotentialMasterCellFasterer(getSpeciesManager(), box, 1, pmBonding.getBondingInfo());
        } else {
            potentialMasterPair = new PotentialMasterFasterer(getSpeciesManager(), box, pmBonding.getBondingInfo());
        }
        potentialMaster = new PotentialComputeAggregate(pmBonding, potentialMasterPair);

        integratorMC = new IntegratorMCFasterer(potentialMaster, random, temperature, box);
        MCMoveRotateArmFasterer rotateArmMove = new MCMoveRotateArmFasterer(potentialMaster, random, box, 1, l);
        integratorMC.getMoveManager().addMCMove(rotateArmMove);
//            ((MCMoveStepTracker) rotateArmMove.getTracker()).setNoisyAdjustment(true);
        MCMoveBondLengthFasterer bondMove = new MCMoveBondLengthFasterer(potentialMaster, random, box, 1, l);
        integratorMC.getMoveManager().addMCMove(bondMove);
//            ((MCMoveStepTracker) bondMove.getTracker()).setNoisyAdjustment(true);

        this.getController().addActivity(new ActivityIntegrate(integratorMC));

        potentialFene = new P2Fene(space);
        potentialWCA = new P2WCA(space);

        P2SoftSphericalSum pBonding = new P2SoftSphericalSum(space, potentialFene, potentialWCA);

        pmBonding.setBondingPotentialPair(species, pBonding, bondedPairs);

        AtomType type = species.getAtomType(0);
        potentialMasterPair.setPairPotential(type, type, potentialWCA);
    }

    public static void main(String[] args) {
        PolymerParam params = new PolymerParam();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        } else {
            params.numSteps = 100000L;
            params.fv = 7;
            params.lv = 34;
        }
        long steps = params.numSteps;
        int f = params.fv;
        int l = params.lv;

        double temperature = 1.0;
        final int dataInterval = 100;
        boolean graphics = true;
        boolean writeConf = false;

        final StarPolymerMCFasterer sim = new StarPolymerMCFasterer(f, l, temperature);

        sim.integratorMC.reset();
        final MeterPotentialEnergyFromIntegratorFasterer meterPE = new MeterPotentialEnergyFromIntegratorFasterer(sim.integratorMC);
        double u = meterPE.getDataAsScalar();
        System.out.println("Initial Potential energy: " + u);

        DataFork forkPE = new DataFork();
        AccumulatorAverageCollapsing accPE = new AccumulatorAverageCollapsing();
        AccumulatorHistory histPE = new AccumulatorHistory(new HistoryCollapsingAverage());
        forkPE.addDataSink(accPE);
        forkPE.addDataSink(histPE);
        DataPumpListener pumpPF = new DataPumpListener(meterPE, forkPE, dataInterval);
        sim.integratorMC.getEventManager().addListener(pumpPF);
        sim.integratorMC.reset();

        MeterRadiusGyration meterRG = new MeterRadiusGyration(sim.getSpace());
        meterRG.setBox(sim.getBox(0));
        System.out.println(String.format("Squared Radius of Gyration: %.8f\n", meterRG.getDataAsScalar()));

        MeterBondLength meterBL = new MeterBondLength(sim.box, f, l);

        if (graphics) {
            final String APP_NAME = "StarPolymer-MC";
//            sim.getBox(0).getBoundary().setBoxSize(Vector.of(new double[]{50, 50, 50}));
            final SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, APP_NAME, 100);
            DisplayBox displayBox = simGraphic.getDisplayBox(sim.getBox(0));
            displayBox.setShowBoundary(false);
            simGraphic.getController().getReinitButton().setPostAction(simGraphic.getPaintAction(sim.box));
            simGraphic.makeAndDisplayFrame(APP_NAME);
            DiameterHashByType dh = (DiameterHashByType) simGraphic.getDisplayBox(sim.box).getDiameterHash();
            dh.setDiameter(sim.species.getAtomType(0), 1.0);
//            sim.ai.setMaxSteps(1);
            DataFork rgFork = new DataFork();
            AccumulatorAverageCollapsing accRG = new AccumulatorAverageCollapsing();

            AccumulatorHistory histRG = new AccumulatorHistory(new HistoryCollapsingAverage());
            rgFork.addDataSink(accRG);
            rgFork.addDataSink(histRG);
            DataPumpListener rgPump = new DataPumpListener(meterRG, rgFork, dataInterval);
            sim.integratorMC.getEventManager().addListener(rgPump);
            sim.integratorMC.reset();

            DisplayPlot pePlot = new DisplayPlot();
            histPE.addDataSink(pePlot.getDataSet().makeDataSink());
            pePlot.setLabel("PE");
            simGraphic.add(pePlot);

            sim.integratorMC.getEventManager().addListener(new IntegratorListenerAction(meterBL, dataInterval));
            DisplayPlot blPlot = new DisplayPlot();
            blPlot.getPlot().setYLog(true);
            DataPumpListener blPump = new DataPumpListener(meterBL, blPlot.getDataSet().makeDataSink(), dataInterval);
            sim.integratorMC.getEventManager().addListener(blPump);
            blPlot.setLabel("bond length");
            simGraphic.add(blPlot);

            DisplayPlot rgPlot = new DisplayPlot();
            histRG.addDataSink(rgPlot.getDataSet().makeDataSink());
            rgPlot.setLabel("Radius of Gyration^2");
            simGraphic.add(rgPlot);
            simGraphic.getController().getDataStreamPumps().add(pumpPF);
            simGraphic.getController().getDataStreamPumps().add(rgPump);
            simGraphic.getController().getResetAveragesButton().setPostAction(new IAction() {
                @Override
                public void actionPerformed() {
                    meterBL.reset();
                }
            });

            return;
        }

        if (writeConf) {
            WriteConfiguration writeConfigurationd = new WriteConfiguration(sim.getSpace());
            writeConfigurationd.setBox(sim.box);
            writeConfigurationd.setConfName("initial conf");
            File file = new File("./resource/init_f" + f);
            writeConfigurationd.setFileName("./resource/init_f" + f);
            writeConfigurationd.actionPerformed();
            sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integratorMC, steps));

            writeConfigurationd.setConfName("new conf");
            File file2 = new File("./resource/newConf_f" + f);
            writeConfigurationd.setFileName("./resource/newConf_f" + f);
            writeConfigurationd.actionPerformed();
            System.out.println("Finished writing configurations");
            System.exit(0);
        }

        System.out.println("Equilibration of " + steps + " steps starts...");
        long t1 = System.currentTimeMillis();
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integratorMC, steps));
System.out.println("Equilibration finished! ");

        DataFork rgFork = new DataFork();
        AccumulatorAverageFixed accRG = new AccumulatorAverageFixed(1000);
        rgFork.addDataSink(accRG);
        DataPumpListener rgPump = new DataPumpListener(meterRG, rgFork, 10);
        sim.integratorMC.getEventManager().addListener(rgPump);
        sim.integratorMC.reset();

        System.out.println("Production starts...");

        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integratorMC, 1000 * 10 * 100));
        System.out.println("Production finished! ");
        long t2 = System.currentTimeMillis();
        System.out.println("time : " + (t2 - t1) / 1000.0);

        System.out.println();
        System.out.println("Final Potential energy: " + meterPE.getDataAsScalar());
        System.out.println(String.format("Final rg^2 : %.8f\n", meterRG.getDataAsScalar()));

        DataGroup rgBase = (DataGroup) accRG.getData();
        IData averageRG = rgBase.getData(accRG.AVERAGE.index);
        IData errorRG = rgBase.getData(accRG.ERROR.index);
        IData correlationRG = rgBase.getData(accRG.BLOCK_CORRELATION.index);
        System.out.print(String.format("Squared RG average: %20.15e error: %9.4e cor: %6.4f\n",
                averageRG.getValue(0), errorRG.getValue(0), correlationRG.getValue(0)));

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
    }

}
