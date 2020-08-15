package etomica.starpolymer;


import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.DiameterHashByType;
import etomica.atom.IAtom;
import etomica.atom.iterator.ApiIndexList;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.history.HistoryCollapsingAverage;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.meter.MeterRadiusGyration;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayPlot;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.*;
import etomica.nbr.CriterionInterMolecular;
import etomica.nbr.NeighborCriterion;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.P2Fene;
import etomica.potential.P2WCA;
import etomica.potential.PotentialGroup;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularNonperiodic;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.species.ISpecies;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

import java.util.ArrayList;

public class StarPolymerMD extends Simulation {

    public Box box;
    public int f, l;
    public SpeciesPolymerMono species;
    public IntegratorBox integrator;
    public IntegratorVelocityVerlet integratorMD;
    public IntegratorMC integratorMC;
    public P2Fene potentialFene;
    public P2WCA potentialWCA;
    public PotentialMaster potentialMaster;

    public StarPolymerMD(int f, int l, double temperature, double tStep, boolean useNbrs) {
        super(Space3D.getInstance());
        this.f = f;
        this.l = l;
        species = new SpeciesPolymerMono(this, getSpace(), f, l);
        species.setIsDynamic(true);
        addSpecies(species);

        box = this.makeBox(new BoundaryRectangularNonperiodic(space));
        box.getBoundary().setBoxSize(new Vector3D(1.5 * l * 2, 1.5 * l * 2, 1.5 * l * 2));
        box.setNMolecules(species, 1);
        if (false) {
            ConformationStarPolymerAll conf = new ConformationStarPolymerAll(this.getSpace(), "./resource/f5L40.xyz", 0);
            conf.initializePositions(box.getMoleculeList().get(0).getChildList());
        }
        ArrayList<ArrayList<int[]>> pairArray = getPairArray(f, l);
        int[][] boundedPairs = pairArray.get(0).toArray(new int[0][0]);
        int[][] nonBoundedPairs = pairArray.get(1).toArray(new int[0][0]);
        boolean doMC = false;
        if (doMC) useNbrs = false;
        if (l < 30) useNbrs = false;
        if (useNbrs) {
            potentialMaster = new PotentialMasterList(this, 2, space);
        } else {
            potentialMaster = new PotentialMaster();
        }
        integratorMD = new IntegratorVelocityVerlet(this, potentialMaster, box);
        integratorMD.setTimeStep(tStep);
        integratorMD.setTemperature(temperature);
        integratorMD.setIsothermal(true);
        integratorMD.setThermostat(IntegratorMD.ThermostatType.ANDERSEN);
        integratorMD.setThermostatInterval(500);

        integratorMC = new IntegratorMC(this, potentialMaster, box);
        MCMoveUpdateConformation moveConformation = new MCMoveUpdateConformation(this, space, "./resource/f5L40.xyz");
        integratorMC.getMoveManager().addMCMove(moveConformation);

        integrator = doMC ? integratorMC : integratorMD;
        this.getController().addActivity(new ActivityIntegrate(integrator));

        potentialFene = new P2Fene(space);
        potentialWCA = new P2WCA(space);

        PotentialGroup potentialGroup = potentialMaster.makePotentialGroup(1);
        ApiIndexList boundedIndexList = new ApiIndexList(boundedPairs);
        ApiIndexList nonBoundedIndexList = new ApiIndexList(nonBoundedPairs);

        potentialGroup.addPotential(potentialFene, boundedIndexList);
        AtomType type = species.getAtomType(0);
        potentialMaster.addPotential(potentialGroup, new ISpecies[]{species});
        if (useNbrs) {
            potentialMaster.addPotential(potentialWCA, new AtomType[]{type, type});
            CriterionInterMolecular criterion = (CriterionInterMolecular) ((PotentialMasterList) potentialMaster).getCriterion(type, type);
//CriterionSimple criterion2 = (CriterionSimple)((CriterionAdapter)criterion.getWrappedCriterion();
//        System.out.println(criterion2);
            criterion.setIntraMolecularCriterion(new NeighborCriterion() {
                @Override
                public boolean accept(IAtom atom1, IAtom atom2) {
//                    int idx1 = atom1.getIndex();
//                    int idx2 = atom2.getIndex();
//                    if (idx1 > idx2) {
//                        int temp = idx1;
//                        idx1 = idx2;
//                        idx2 = temp;
//                    }
//                    if (idx1 == 0 && idx2 % l == 1) return false;
//                    if (idx2 - 1 == idx1 && idx1 % l != 0) return false;
                    return true;
                }

                @Override
                public boolean needUpdate(IAtom atom) {
                    return criterion.needUpdate(atom);
                }

                @Override
                public void setBox(Box box) {
                    criterion.setBox(box);
                }

                @Override
                public boolean unsafe() {
                    return criterion.unsafe();
                }

                @Override
                public void reset(IAtom atom) {
                    criterion.reset(atom);
                }
            });
            criterion.setIntraMolecularOnly(true);
            integrator.getEventManager().addListener(((PotentialMasterList) potentialMaster).getNeighborManager(box));
        } else {
            potentialGroup.addPotential(potentialWCA, nonBoundedIndexList);
            potentialGroup.addPotential(potentialWCA, boundedIndexList);
        }
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

        final MeterPotentialEnergy meterPE = new MeterPotentialEnergy(sim.potentialMaster, sim.box);
        double u = meterPE.getDataAsScalar();
        System.out.println("Initial Potential energy: " + u);

        DataFork forkPE = new DataFork();
        AccumulatorAverageCollapsing accPE = new AccumulatorAverageCollapsing();
        AccumulatorHistory histPE = new AccumulatorHistory(new HistoryCollapsingAverage());
        forkPE.addDataSink(accPE);
        forkPE.addDataSink(histPE);
        DataPumpListener pumpPF = new DataPumpListener(meterPE, forkPE, dataInterval);
        sim.integrator.getEventManager().addListener(pumpPF);
        sim.integrator.reset();

        MeterRadiusGyration meterRG = new MeterRadiusGyration(sim.getSpace());
        meterRG.setBox(sim.getBox(0));
        System.out.println(String.format("Squared Radius of Gyration: %.8f\n", meterRG.getDataAsScalar()));

        DataFork rgFork = new DataFork();
        AccumulatorAverageCollapsing accRG = new AccumulatorAverageCollapsing();
        AccumulatorHistory histRG = new AccumulatorHistory(new HistoryCollapsingAverage());
        rgFork.addDataSink(accRG);
        rgFork.addDataSink(histRG);
        DataPumpListener rgPump = new DataPumpListener(meterRG, rgFork, dataInterval);
        sim.integrator.getEventManager().addListener(rgPump);
        sim.integrator.reset();

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

            sim.integrator.getEventManager().addListener(new IntegratorListenerAction(meterBL, dataInterval));
            DisplayPlot blPlot = new DisplayPlot();
            blPlot.getPlot().setYLog(true);
            DataPumpListener blPump = new DataPumpListener(meterBL, blPlot.getDataSet().makeDataSink(), dataInterval);
            sim.integrator.getEventManager().addListener(blPump);
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
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator), steps);
        long t2 = System.currentTimeMillis();
        System.out.println("MD finished! ");
        System.out.println("time : " + (t2 - t1) / 1000.0);

        if (true) {
            System.out.println();
            System.out.println("MD begin for dumping rg^2 ");
            int interval = 1000;
            DataLogger dataLogger = new DataLogger();
            DataPumpListener rgPumpListener = new DataPumpListener(meterRG, dataLogger, interval);
            sim.integrator.getEventManager().addListener(rgPumpListener);
            dataLogger.setFileName(fp);
            dataLogger.setAppending(false);
            DataArrayWriter writer = new DataArrayWriter();
            writer.setIncludeHeader(false);
            dataLogger.setDataSink(writer);
            long t3 = System.currentTimeMillis();
            sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator), xsteps);
            dataLogger.cleanUp();
            System.out.println("Dumping finished! ");
            System.out.println("time: " + (t3 - t1) / 1000.0);
        }
        System.out.println();
        System.out.println("Final Potential energy: " + meterPE.getDataAsScalar());
        System.out.println(String.format("Final rg^2 : %.8f\n", meterRG.getDataAsScalar()));

    }

    private ArrayList<ArrayList<int[]>> getPairArray(int f, int l) {
        ArrayList<ArrayList<int[]>> pairArray = new ArrayList<>(2);
        ArrayList<int[]> boundedPairArray = new ArrayList<>();
        ArrayList<int[]> nonBoundedPariArray = new ArrayList<>();

        for (int k = 1; k < f * l + 1; k++) {
            int[] temp = new int[]{0, k};

            if (k % l == 1) {
                boundedPairArray.add(temp);
            } else {
                nonBoundedPariArray.add(temp);
            }
        }
        for (int i = 1; i < f * l + 1; i++) {
            for (int j = i + 1; j < f * l + 1; j++) {
                int[] temp = new int[]{i, j};
                if (j != i + 1 | i % l == 0) {
                    nonBoundedPariArray.add(temp);
                } else {
                    boundedPairArray.add(temp);
                }
            }
        }
        pairArray.add(boundedPairArray);
        pairArray.add(nonBoundedPariArray);
        return pairArray;
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
