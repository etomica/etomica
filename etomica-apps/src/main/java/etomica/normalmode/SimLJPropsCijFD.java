package etomica.normalmode;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.history.HistoryCollapsingAverage;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.meter.MeterPressure;
import etomica.graphics.ColorScheme;
import etomica.graphics.DisplayPlot;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.lattice.crystal.*;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.*;
import etomica.simulation.Simulation;
import etomica.space.*;
import etomica.species.SpeciesGeneral;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.util.random.RandomMersenneTwister;

import java.awt.*;

public class SimLJPropsCijFD extends Simulation {
    public CoordinateDefinitionLeaf coordinateDefinition, coordinateDefinitionp, coordinateDefinitionm;
    public IntegratorVelocityVerlet integrator;
    public Box box, boxP, boxM;
    public Boundary boundary, boundaryp, boundarym;
    public int[] nCells;
    public Basis basis, basisp, basism;
    public Primitive primitive, primitivep, primitivem;
    public PotentialMasterList potentialMaster, potentialMasterP, potentialMasterM;
    public P2SoftSphericalTruncated potential;
    public SpeciesGeneral species;

    public SimLJPropsCijFD(Space _space, int N, double rho, double T, double rc, double dt, double deV, boolean isLRC) {
        super(_space);

        setRandom(new RandomMersenneTwister(1)); // random seed

        species = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this), true);
        addSpecies(species);
        double L = Math.pow(4.0 / rho, 1.0 / 3.0);
        int nc = (int) Math.round(Math.pow(N / 4, 1.0 / 3.0));
        boundary = new BoundaryRectangularPeriodic(space, nc * L);
        box = this.makeBox(boundary);
        box.setNMolecules(species, N);
        primitive = new PrimitiveCubic(space, nc * L);
        nCells = new int[]{nc, nc, nc};
        Basis basisFCC = new BasisCubicFcc();
        basis = new BasisBigCell(space, basisFCC, nCells);
        coordinateDefinition = new CoordinateDefinitionLeaf(box, primitive, basis, space);
        coordinateDefinition.initializeCoordinates(new int[]{1, 1, 1});
        potentialMaster = new PotentialMasterList(getSpeciesManager(), box, 2, rc, BondingInfo.noBonding());
        potentialMaster.doAllTruncationCorrection = isLRC;

        integrator = new IntegratorVelocityVerlet(potentialMaster, getRandom(), dt, T, box);
        integrator.setIsothermal(true);
        integrator.setThermostatNoDrift(true);
        potential = new P2SoftSphericalTruncated(new P2LennardJones(), rc);
        AtomType sphereType = species.getLeafType();
        potentialMaster.setPairPotential(sphereType, sphereType, potential);
        potentialMaster.init(); // should not need to do this in MeterSolidHMAFD
        integrator.getEventManager().removeListener(potentialMaster);
        this.getController().addActivity(new ActivityIntegrate(integrator));

        if (!isLRC) {
            ((P2SoftSphericalTruncated) potential).setTruncationRadius(0.6 * boundary.getBoxSize().getX(0));
        }

        //box+
        double rhop = rho / (1 + deV);
        double Lp = Math.pow(4.0/rhop, 1.0 / 3.0);
        boundaryp = new BoundaryRectangularPeriodic(space, nc * Lp);
        boxP = this.makeBox(boundaryp);
        boxP.setNMolecules(species, N);
        primitivep = new PrimitiveCubic(space, nc * Lp);
        Basis basisFCCp = new BasisCubicFcc();
        basisp = new BasisBigCell(space, basisFCCp, nCells);
        coordinateDefinitionp = new CoordinateDefinitionLeaf(boxP, primitivep, basisp, space);
        coordinateDefinitionp.initializeCoordinates(new int[]{1, 1, 1});
        potentialMasterP = new PotentialMasterList(getSpeciesManager(), boxP, 2, rc, BondingInfo.noBonding());
        potentialMasterP.doAllTruncationCorrection = false;
        potentialMasterP.setPairPotential(sphereType, sphereType, potential);
        potentialMasterP.init();

        //box-
        double rhom = rho / (1 - deV);
        double Lm = Math.pow(4.0/rhom, 1.0 / 3.0);
        boundarym = new BoundaryRectangularPeriodic(space, nc * Lm);
        boxM = this.makeBox(boundarym);
        boxM.setNMolecules(species, N);
        primitivem = new PrimitiveCubic(space, nc * Lm);
        Basis basisFCCm = new BasisCubicFcc();
        basism = new BasisBigCell(space, basisFCCm, nCells);
        coordinateDefinitionm = new CoordinateDefinitionLeaf(boxM, primitivem, basism, space);
        coordinateDefinitionm.initializeCoordinates(new int[]{1, 1, 1});
        potentialMasterM = new PotentialMasterList(getSpeciesManager(), boxM, 2, rc, BondingInfo.noBonding());
        potentialMasterM.doAllTruncationCorrection = false;
        potentialMasterM.setPairPotential(sphereType, sphereType, potential);
        potentialMasterM.init();

    }

    public void initialize(long initSteps) {
        this.getController().runActivityBlocking(new ActivityIntegrate(this.integrator, initSteps));
    }

    public static void main(String[] args) {
        SimLJPropsCijParam params = new SimLJPropsCijParam();
        ParseArgs.doParseArgs(params, args);
        long steps = params.steps;
        long stepsEq = params.stepsEq;
        long stepsShort = params.stepsEq;
        int N = params.N;
        double rho = params.rho;
        double rc = params.rc;
        double dbeta = params.dbeta;
        double deV = params.deV;
        boolean isLRC = params.isLRC;
        double T = params.T;
        double dt = params.dt;
        double gV = params.gV;
        double gVV = params.gVV;
        boolean isGraphic = params.isGraphic;
        double[] elasticParams = new double[11];
        elasticParams[0] = gV;
        elasticParams[1] = gVV;


        final SimLJPropsCijFD sim = new SimLJPropsCijFD(Space.getInstance(3), N, rho, T, rc, dt, deV, isLRC);
        double V = sim.box.getBoundary().volume();

        //Ulat
        sim.integrator.reset();//uLat
        MeterPotentialEnergyFromIntegrator meterPE = new MeterPotentialEnergyFromIntegrator(sim.integrator);
        double uLat = meterPE.getDataAsScalar();

        //Plat
        MeterPressure meterP = new MeterPressure(sim.box, sim.potentialMaster);
        meterP.setTemperature(0);//ZERO means NO ideal gas component!
        double pLat = meterP.getDataAsScalar();
        meterP.setTemperature(T);


        System.out.println("N " + N + " V " + V + " L " + Math.pow(V, 1.0 / 3.0));
        System.out.println("rho " + rho + " T " + T);
        System.out.println("rc " + rc);
        System.out.println("steps " + steps + " stepsEq " + stepsEq + " stepsShort " + stepsShort);
        System.out.println("dt " + dt);

        System.out.println("\nHMA Elastic Parameters");
        System.out.println("gV  " + params.gV + " gVV  " + params.gVV);

        System.out.println("\nLattice properties");
        System.out.println("uLat  " + uLat);
        System.out.println("pLat  " + pLat);

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
            simGraphic.makeAndDisplayFrame(("LJ-FCC"));

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

            simGraphic.add(historyPE);
            simGraphic.add(historyP);
            return;
        }

        MeterSolidHMAFD meterElastic = new MeterSolidHMAFD(sim.space, new MeterPotentialEnergyFromIntegrator(sim.integrator), sim.potentialMaster,sim.potentialMasterP,sim.potentialMasterM, sim.boxP, sim.boxM, sim.coordinateDefinition, T, elasticParams, dbeta, deV);

//Equilibration
        final long startTime = System.currentTimeMillis();
        System.out.flush();
        sim.initialize(stepsEq);

        int numBlocks = 100;
        int interval = 10;
        long blockSize = steps / (numBlocks * interval);
        if (blockSize == 0) blockSize = 1;
        System.out.println("block size " + blockSize + "  interval " + interval);


        //Mapped
        AccumulatorAverageCovariance accumulatorElastic = new AccumulatorAverageCovariance(blockSize);
        DataPumpListener accumulatorUPPump = new DataPumpListener(meterElastic, accumulatorElastic, interval);
        sim.integrator.getEventManager().addListener(accumulatorUPPump);

        System.out.println("Short sim for Covariance  " + stepsShort + " steps");
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, stepsShort));

        IData dataElasticAvg0 = accumulatorElastic.getData(accumulatorElastic.AVERAGE);
        IData dataElasticErr0 = accumulatorElastic.getData(accumulatorElastic.ERROR);
        double uShift = dataElasticAvg0.getValue(1);
        double errUshift = dataElasticErr0.getValue(1);
        double pShift = dataElasticAvg0.getValue(3);
        double errPshift = dataElasticErr0.getValue(3);

        accumulatorElastic.reset();
        System.out.println("uShift  " + uShift + "  " + errUshift);
        System.out.println("pShift  " + pShift + "  " + errPshift);

        meterElastic.setShift(uShift, pShift);

// Production run ...
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, steps));

        IData dataElasticAvg = accumulatorElastic.getData(accumulatorElastic.AVERAGE);
        IData dataElasticErr = accumulatorElastic.getData(accumulatorElastic.ERROR);
        IData dataElasticCor = accumulatorElastic.getData(accumulatorElastic.BLOCK_CORRELATION);
        IData dataElasticCov = accumulatorElastic.getData(accumulatorElastic.COVARIANCE);


/** First Derivatives*/
        double Uc     = dataElasticAvg.getValue(0) + uShift;
        double errUc  = dataElasticErr.getValue(0);
        double corUc = dataElasticCor.getValue(0);

        double Um      = dataElasticAvg.getValue(1) + uShift;
        double errUm   = dataElasticErr.getValue(1);
        double corUm  = dataElasticCor.getValue(1);

        double Pc     = dataElasticAvg.getValue(2) + pShift;
        double errPc  = dataElasticErr.getValue(2);
        double corPc = dataElasticCor.getValue(2);

        double Pm      = dataElasticAvg.getValue(3) + pShift;
        double errPm   = dataElasticErr.getValue(3);
        double corPm  = dataElasticCor.getValue(3);


        /** Second Derivatives*/
        int nd = dataElasticAvg.getLength();
        double varUc  = dataElasticCov.getValue(0 * (nd + 1)) * N * N;
        double varUm   = dataElasticCov.getValue(1 * (nd + 1)) * N * N;
        double varPc  = dataElasticCov.getValue(2 * (nd + 1));
        double varPm   = dataElasticCov.getValue(3 * (nd + 1));
        double covUcPc  = dataElasticCov.getValue(0 * nd + 2) * N;
        double covUmPm   = dataElasticCov.getValue(1 * nd + 3) * N;

        double T2 = T * T;

        double Cvc = varUc / T2 / N;
        double Cvm  = varUm / T2 / N - dataElasticAvg.getValue(4) + 3.0 / 2.0 * (N - 1.0) / N ;

        double Bc = dataElasticAvg.getValue(5) - V / T * varPc;
        double Bm = dataElasticAvg.getValue(6) - V / T * varPm;

        double Gc = rho + covUcPc / T2;
        double Gm = dataElasticAvg.getValue(7) + covUmPm / T2 + Pm / T;


        //Print
        System.out.println("\n++++ First Derivatives ++++");
        System.out.println("Uc   " + Uc + " +/- " + errUc + "   cor " + corUc);
        System.out.println("Um    " + Um + " +/- " + errUm + "   cor " + corUm);
        System.out.println("Pc   " + Pc + " +/- " + errPc + "   cor " + corPc);
        System.out.println("Pm    " + Pm + " +/- " + errPm + "   cor " + corPm);

        System.out.println("\n ++++ Second Derivatives ++++");
        System.out.println("Cvc " + Cvc);
        System.out.println("Cvc " + Cvm);
        System.out.println("Bc  " + Bc);
        System.out.println("Bm  " + Bm);
        System.out.println("Gc  " + Gc);
        System.out.println("Gm  " + Gm);

        long endTime = System.currentTimeMillis();
        System.out.println("\n Time " + (endTime - startTime) / 1000.0);
    }

    public static class SimLJPropsCijParam extends ParameterBase {
        public int N = 500;
        public double rho = 1.0;
        public long steps = 10000;
        public long stepsEq = steps / 10;
        public double T = 1;
        public double rc = 3;
        public double dt = 0.002;
        public double gV = 3.094707617816900;
        public double gVV = -0.199514733065326;
        public double dbeta = 0.001;
        public double deV = 0.001;
        public boolean isGraphic = !true;
        public boolean isLRC = false;
    }
}