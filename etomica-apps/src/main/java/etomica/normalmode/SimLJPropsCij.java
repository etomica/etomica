package etomica.normalmode;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.history.HistoryCollapsingAverage;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.meter.MeterPressure;
import etomica.data.meter.MeterPressureTensor;
import etomica.data.types.DataGroup;
import etomica.graphics.ColorScheme;
import etomica.graphics.DisplayPlot;
import etomica.graphics.DisplayTextBox;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisCubicFcc;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.potential.Potential2SoftSpherical;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.species.SpeciesGeneral;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

import java.awt.*;

public class SimLJPropsCij extends Simulation {
    public final CoordinateDefinitionLeaf coordinateDefinition;
    public IntegratorMC integrator;
    public ActivityIntegrate activityIntegrate;
    public Box box;
    public Boundary boundary;
    public int[] nCells;
    public Basis basis;
    public Primitive primitive;
    public MCMoveAtomCoupled atomMove;
    public PotentialMasterList potentialMaster;
    public Potential2SoftSpherical potential;
    public SpeciesGeneral species;

    public SimLJPropsCij(Space _space, int numAtoms, double density, double temperature, double rc, boolean isLRC, double ex, double gamma) {
        super(_space);

        species = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this));
        addSpecies(species);


        potentialMaster = new PotentialMasterList(this, space);

        // TARGET
        double L = Math.pow(4.0 / density, 1.0 / 3.0);
        int n = (int) Math.round(Math.pow(numAtoms / 4, 1.0 / 3.0));
        boundary = new BoundaryRectangularPeriodic(space, n * L);
        box = this.makeBox(boundary);
        box.setNMolecules(species, numAtoms);

        integrator = new IntegratorMC(potentialMaster, getRandom(), temperature, box);
        MeterPotentialEnergy meterPE = new MeterPotentialEnergy(potentialMaster, box);
        atomMove = new MCMoveAtomCoupled(potentialMaster, meterPE, getRandom(), space);
        atomMove.setStepSize(0.1);
        atomMove.setStepSizeMax(0.5);
        atomMove.setDoExcludeNonNeighbors(true);
        integrator.getMoveManager().addMCMove(atomMove);
//        ((MCMoveStepTracker)atomMove.getTracker()).setNoisyAdjustment(true);

        primitive = new PrimitiveCubic(space, n * L);

        nCells = new int[]{n, n, n};
        Basis basisFCC = new BasisCubicFcc();
        basis = new BasisBigCell(space, basisFCC, nCells);

        coordinateDefinition = new CoordinateDefinitionLeaf(box, primitive, basis, space);
        coordinateDefinition.initializeCoordinates(new int[]{1, 1, 1});

        Potential2SoftSpherical potential = new P2LennardJones(space, 1.0, 1.0);
        potential = new P2SoftSphericalTruncated(space, potential, rc);
        atomMove.setPotential(potential);

        AtomType sphereType = species.getLeafType();
        potentialMaster.addPotential(potential, new AtomType[]{sphereType, sphereType});

        potentialMaster.lrcMaster().setEnabled(false);

        int cellRange = 2;
        potentialMaster.setRange(rc);
        potentialMaster.setCellRange(cellRange); // NeighborCellManager handles this even if cells are a bit small
        // find neighbors now.  Don't hook up NeighborListManager (neighbors won't change)
        potentialMaster.getNeighborManager(box).reset();
        int potentialCells = potentialMaster.getNbrCellManager(box).getLattice().getSize()[0];
        if (potentialCells < cellRange * 2 + 1) {
            throw new RuntimeException("oops (" + potentialCells + " < " + (cellRange * 2 + 1) + ")");
        }

        getController().addActivity(new ActivityIntegrate(integrator));

        ((P2SoftSphericalTruncated) potential).setTruncationRadius(0.6 * boundary.getBoxSize().getX(0));

    }
    // End of Constructor

    public void initialize(long initSteps) {
        this.getController().runActivityBlocking(new ActivityIntegrate(this.integrator, initSteps));
        integrator.getMoveManager().setEquilibrating(false);
    }


    public static void main(String[] args) {
        SimOverlapParam params = new SimOverlapParam();
        ParseArgs.doParseArgs(params, args);
        long numSteps = params.numSteps;
        final int numAtoms = params.numAtoms;
        double density = params.density;
        double rc = params.rc;
        double temperature = params.temperature;
        double dbP = params.dbP;
        double dU_est = params.dU_est;
        double ddbP = params.ddbP;
        double dP = dbP * temperature;
        double ddP = -(ddbP * density * density * temperature) / numAtoms; // ddP is dP/dV NOT dP/drho
        boolean isLRC = params.isLRC;
        boolean isGraphic = params.isGraphic;
        double dSigma11 = params.dSigma11;
        double dc11 = params.dc11;
        double gamma = params.gamma;
        double ex = params.ex;

        System.out.println(numAtoms + " atoms at density " + density + " and temperature " + temperature + " rc = " + rc);
        System.out.println(numSteps + " MC steps");
        System.out.println(" ex = " + ex + " gamma = " + gamma + " dSigma11 = " + dSigma11 + " dc11 = " + dc11);

        final SimLJPropsCij sim = new SimLJPropsCij(Space.getInstance(3), numAtoms, density, temperature, rc, isLRC, ex, gamma);



        MeterPotentialEnergyFromIntegrator meterPE = new MeterPotentialEnergyFromIntegrator(sim.integrator);
        double ULat = meterPE.getDataAsScalar();
        //P
        MeterPressure meterP = new MeterPressure(sim.space);
        meterP.setBox(sim.box);
        meterP.setTemperature(0);//ZERO means NO ideal gas component!
        meterP.setPotentialMaster(sim.potentialMaster);
        double PLat = meterP.getDataAsScalar();
        meterP.setTemperature(temperature);

        //Sigma
        MeterPressureTensor meterPij = new MeterPressureTensor(sim.potentialMaster, sim.space);
        meterPij.setBox(sim.box);
        meterPij.setTemperature(0);
//        meterSigma.getPotentialCalculation().setIntegrator(sim.integrator);
        IData dataSigma11_lat = meterPij.getData();
        meterPij.setTemperature(temperature);
        dataSigma11_lat.TE(-1);
        double sigma11_lat = dataSigma11_lat.getValue(0);
        System.out.println("sigma11_lat = " + sigma11_lat);

        System.out.println(" Lattice Energy/N = " + (ULat / numAtoms));
        System.out.println(" Lattice Pressure = " + PLat);


        System.out.println(" Uharm = " + 1.5 * (numAtoms - 1.0) / numAtoms * temperature + "  Pharm = " + dP + " dPharm/dRho = " + ddP);// This is dP/dV NOT dP/drho
//        System.out.println("dU_est = " + dU_est);
        double volume = sim.box.getBoundary().volume();
        System.out.println(" volume = " + volume);
        MeterSolidPropsLJ meterUP = new MeterSolidPropsLJ(sim.space, new MeterPotentialEnergyFromIntegrator(sim.integrator),
                sim.potentialMaster, sim.coordinateDefinition, temperature, dP, dU_est, ddP, ULat, PLat, dSigma11, dc11, gamma, ex);












        if (isGraphic) {
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
            simGraphic.setPaintInterval(sim.box, 1000);
            ColorScheme colorScheme = new ColorScheme() {
                protected Color[] allColors;
                public Color getAtomColor(IAtom a) {
                    if (allColors == null) {
                        allColors = new Color[768];
                        for (int i = 0; i < 256; i++) {
                            allColors[i] = new Color(255 - i, i, 0);
                        }
                        for (int i = 0; i < 256; i++) {
                            allColors[i + 256] = new Color(0, 255 - i, i);
                        }
                        for (int i = 0; i < 256; i++) {
                            allColors[i + 512] = new Color(i, 0, 255 - i);
                        }
                    }
                    return allColors[(2 * a.getLeafIndex()) % 768];
                }
            };
            simGraphic.getDisplayBox(sim.box).setColorScheme(colorScheme);
            simGraphic.makeAndDisplayFrame(("LJ") + " FCC");



            //PE
            DataSourceCountSteps timeSource = new DataSourceCountSteps(sim.integrator);
            AccumulatorHistory accPE = new AccumulatorHistory(new HistoryCollapsingAverage());
            accPE.setTimeDataSource(timeSource);
            DataPumpListener pumpPE = new DataPumpListener(meterPE, accPE, 10);
            sim.integrator.getEventManager().addListener(pumpPE);
            DisplayPlot historyPE = new DisplayPlot();
            accPE.setDataSink(historyPE.getDataSet().makeDataSink());
            historyPE.setLabel("PE");

            //Pressure
            AccumulatorHistory accP = new AccumulatorHistory(new HistoryCollapsingAverage());
            accP.setTimeDataSource(timeSource);
            DataPumpListener pumpP = new DataPumpListener(meterP, accP, 10);
            sim.integrator.getEventManager().addListener(pumpP);
            DisplayPlot historyP = new DisplayPlot();
            accP.setDataSink(historyP.getDataSet().makeDataSink());
            historyP.setLabel("P");


            // P11
            DisplayPlot historyPij = new DisplayPlot();
            historyPij.setLabel("Pij");

            // Conv P11
            AccumulatorHistory accPij = new AccumulatorHistory(new HistoryCollapsingAverage());
            accPij.setTimeDataSource(timeSource);
            DataSplitter splitter = new DataSplitter();
            DataPumpListener pumpPij = new DataPumpListener(meterPij, splitter, 10);
            splitter.setDataSink(0, accPij);
            sim.integrator.getEventManager().addListener(pumpPij);
            accPij.setDataSink(historyPij.getDataSet().makeDataSink());
            historyPij.setLegend(new DataTag[]{accPij.getTag()}, "Conv");

            // Mapped P11
            DataSplitter splitterM = new DataSplitter();
            DataPumpListener pumpPijM = new DataPumpListener(meterUP, splitterM, 10);
            //along dr
            AccumulatorHistory accPijM1 = new AccumulatorHistory(new HistoryCollapsingAverage());
            accPijM1.setTimeDataSource(timeSource);
            splitterM.setDataSink(11, accPijM1);
            sim.integrator.getEventManager().addListener(pumpPijM);
            accPijM1.setDataSink(historyPij.getDataSet().makeDataSink());
            historyPij.setLegend(new DataTag[]{accPijM1.getTag()}, "HMA: along dr  ");

            //along dr0
            AccumulatorHistory accPijM2 = new AccumulatorHistory(new HistoryCollapsingAverage());
            accPijM2.setTimeDataSource(timeSource);
            splitterM.setDataSink(12, accPijM2);
            accPijM2.setDataSink(historyPij.getDataSet().makeDataSink());
            historyPij.setLegend(new DataTag[]{accPijM2.getTag()}, "HMA: along dr0");


            simGraphic.add(historyPE);
            simGraphic.add(historyP);
            simGraphic.add(historyPij);

            return;
        }












//Initialization
        System.out.flush();
        long Ninit = numSteps / 5;
        sim.initialize(Ninit);

        int numBlocks = 100;
        int interval = numAtoms;
        long blockSize = numSteps / (numBlocks * interval);
        if (blockSize == 0) blockSize = 1;
        System.out.println("block size " + blockSize + " interval " + interval + "\n");

        //U
        AccumulatorAverageFixed accumulatorPE = new AccumulatorAverageFixed(blockSize);
        DataPumpListener accumulatorPEPump = new DataPumpListener(meterPE, accumulatorPE, interval);
        sim.integrator.getEventManager().addListener(accumulatorPEPump);
        //P
        AccumulatorAverageFixed accumulatorP = new AccumulatorAverageFixed(blockSize);
        DataPumpListener accumulatorPPump = new DataPumpListener(meterP, accumulatorP, interval);
        sim.integrator.getEventManager().addListener(accumulatorPPump);

        //Sigma
        AccumulatorAverageFixed accumulatorSigma = new AccumulatorAverageFixed(blockSize);
        DataPumpListener accumulatorSigmaPump = new DataPumpListener(meterPij, accumulatorSigma, interval);
        sim.integrator.getEventManager().addListener(accumulatorSigmaPump);

        //Mapped
        AccumulatorAverageFixed accumulatorUP = new AccumulatorAverageFixed(blockSize);
        DataPumpListener accumulatorUPPump = new DataPumpListener(meterUP, accumulatorUP, interval);
        sim.integrator.getEventManager().addListener(accumulatorUPPump);


        final long startTime = System.currentTimeMillis();
/** RUN ... */
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, numSteps));


        //U

        DataGroup dataPE = (DataGroup) accumulatorPE.getData();
        IData dataPEAvg = dataPE.getData(accumulatorPE.AVERAGE.index);
        IData dataPEErr = dataPE.getData(accumulatorPE.ERROR.index);
        IData dataPECorrelation = dataPE.getData(accumulatorPE.BLOCK_CORRELATION.index);
        double peAvg = dataPEAvg.getValue(0);
        double peErr = dataPEErr.getValue(0);
        double peCor = dataPECorrelation.getValue(0);
        System.out.println(" Udirect/N = " + peAvg / numAtoms + "  Err = " + peErr / numAtoms + "  cor: " + peCor);
        //P
        DataGroup dataP = (DataGroup) accumulatorP.getData();
        IData dataPAvg = dataP.getData(accumulatorP.AVERAGE.index);
        IData dataPErr = dataP.getData(accumulatorP.ERROR.index);
        IData dataPCorrelation = dataP.getData(accumulatorP.BLOCK_CORRELATION.index);
        double pAvg = dataPAvg.getValue(0);
        double pErr = dataPErr.getValue(0);
        double pCor = dataPCorrelation.getValue(0);
        System.out.println(" Pdirect = " + pAvg + "  Err = " + pErr + "  cor: " + pCor);


        //Sigma
        DataGroup dataSigma11 = (DataGroup) accumulatorSigma.getData();
        IData dataSigma11Avg = dataSigma11.getData(accumulatorSigma.AVERAGE.index);
        IData dataSigma11Err = dataSigma11.getData(accumulatorSigma.ERROR.index);
        IData dataSigma11Correlation = dataSigma11.getData(accumulatorSigma.BLOCK_CORRELATION.index);

        density = numAtoms / volume;

        double sigma11Avg = dataSigma11Avg.getValue(0);
        double sigma11Err = dataSigma11Err.getValue(0);
        double sigma11Cor = dataSigma11Correlation.getValue(0);
        System.out.println(" Sigma11Direct= " + sigma11Avg + "  Err = " + sigma11Err + "  cor: " + sigma11Cor);


        DataGroup dataUP = (DataGroup) accumulatorUP.getData();
        IData dataUPAvg = dataUP.getData(accumulatorUP.AVERAGE.index);
        IData dataUPErr = dataUP.getData(accumulatorUP.ERROR.index);
        IData dataUPCorr = dataUP.getData(accumulatorUP.BLOCK_CORRELATION.index);

        double dU = dataUPAvg.getValue(0);
        double Ur = dataUPAvg.getValue(1); //U + 0.5*fdr
        double Pvir = dataUPAvg.getValue(2);   //fr/3/volume
        double Pr = dataUPAvg.getValue(3);
        System.out.println();
        System.out.println("************************************************************************");


        System.out.println(" U/N   " + (ULat + dU) / numAtoms + "  +/- " + dataUPErr.getValue(0) / numAtoms + "  corr: " + dataUPCorr.getValue(0) + " dU = " + dU / numAtoms);
        System.out.println(" Um/N  " + (ULat + 1.5 * (numAtoms - 1) * temperature + Ur) / numAtoms + "  +/- "
                + dataUPErr.getValue(1) / numAtoms + "  corr: " + dataUPCorr.getValue(1) + "   Uanh = " + Ur / numAtoms);

        System.out.println("************************************************************************");
        System.out.println(" P  " + (density * temperature + Pvir) + "  +/- " + dataUPErr.getValue(2) + "  corr: " + dataUPCorr.getValue(2));
        System.out.println(" Pm " + (dP + Pr) + "  +/- " + dataUPErr.getValue(3) + "  corr: " + dataUPCorr.getValue(3));
        System.out.println("************************************************************************");


        double sigma11_conv =  dataUPAvg.getValue(10);
        double sigma11r  = dataUPAvg.getValue(11);
        double sigma11r2 = dataUPAvg.getValue(12);

        System.out.println(" sigma11_conv    : " + sigma11_conv + " +/- " + dataUPErr.getValue(10) + "   corr: " + dataUPCorr.getValue(10) +
                "    dSigma11_conv: " + (sigma11_conv - sigma11_lat));
        //Mapping along dr
        System.out.println(" sigma11_mapped1 : " + sigma11r + " +/- " + dataUPErr.getValue(11) + "   corr: " + dataUPCorr.getValue(11));
        //Mapping along dr0
        System.out.println(" sigma11_mapped2 : " + sigma11r2 + " +/- " + dataUPErr.getValue(12) + "   corr: " + dataUPCorr.getValue(12));
        //Difference from Conv
        System.out.println("      diff1 = " + (sigma11r - sigma11_conv) + "   diff2 = " + (sigma11r2 - sigma11_conv));

        long endTime = System.currentTimeMillis();
        System.out.println("time: " + (endTime - startTime) / 1000.0);
    }


    public static class SimOverlapParam extends ParameterBase {
        int n = 5;
        public int numAtoms = 4 * n * n * n;
        public double density = 1.0;
        public long numSteps = 1000000;
        public double temperature = 1.0; //Tm=0.930 at rho=1.0
        public double dbP = 9.6 / temperature;
        public boolean isLRC = false;
        public boolean isGraphic = true;
        public double dU_est = 272.814649;
        public double ddbP = 6.933;
        public double rc = 2.5;
        public double ex = 0.0;
        public double gamma = 0.0;
        public double dSigma11 = 1.4;
        public double dc11 = 0.7542631471006018;
    }
}