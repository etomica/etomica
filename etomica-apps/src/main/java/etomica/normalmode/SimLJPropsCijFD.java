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

import java.awt.*;

public class SimLJPropsCijFD extends Simulation {
    public CoordinateDefinitionLeaf coordinateDefinition, coordinateDefinitionp, coordinateDefinitionm;
    public IntegratorVelocityVerlet integrator;
    public Box box, boxp, boxm;
    public Boundary boundary, boundaryp, boundarym;
    public int[] nCells;
    public Basis basis, basisp, basism;
    public Primitive primitive, primitivep, primitivem;
    public PotentialMasterList potentialMaster, potentialMasterp, potentialMasterm;
    public P2SoftSphericalTruncated potential;
    public SpeciesGeneral species;

    public SimLJPropsCijFD(Space _space, int N, double rho, double T, double rc, double dt, double deV, boolean isLRC) {
        super(_space);
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
        boxp = this.makeBox(boundaryp);
        boxp.setNMolecules(species, N);
        primitivep = new PrimitiveCubic(space, nc * Lp);
        Basis basisFCCp = new BasisCubicFcc();
        basisp = new BasisBigCell(space, basisFCCp, nCells);
        coordinateDefinitionp = new CoordinateDefinitionLeaf(boxp, primitivep, basisp, space);
        coordinateDefinitionp.initializeCoordinates(new int[]{1, 1, 1});
        potentialMasterp = new PotentialMasterList(getSpeciesManager(), boxp, 2, rc, BondingInfo.noBonding());
        potentialMasterp.doAllTruncationCorrection = false;
        potentialMasterp.setPairPotential(sphereType, sphereType, potential);
        potentialMasterp.init();

        //box-
        double rhom = rho / (1 - deV);
        double Lm = Math.pow(4.0/rhom, 1.0 / 3.0);
        boundarym = new BoundaryRectangularPeriodic(space, nc * Lm);
        boxm = this.makeBox(boundarym);
        boxm.setNMolecules(species, N);
        primitivem = new PrimitiveCubic(space, nc * Lm);
        Basis basisFCCm = new BasisCubicFcc();
        basism = new BasisBigCell(space, basisFCCm, nCells);
        coordinateDefinitionm = new CoordinateDefinitionLeaf(boxm, primitivem, basism, space);
        coordinateDefinitionm.initializeCoordinates(new int[]{1, 1, 1});
        potentialMasterm = new PotentialMasterList(getSpeciesManager(), boxm, 2, rc, BondingInfo.noBonding());
        potentialMasterm.doAllTruncationCorrection = false;
        potentialMasterm.setPairPotential(sphereType, sphereType, potential);
        potentialMasterm.init();

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


        System.out.println(" N " + N + " volume " + V);
        System.out.println(" T " + T + " rho " + rho );
        System.out.println(" rc " + rc);
        System.out.println(" steps " + steps + " stepsEq " + stepsEq + " stepsShort = " + stepsShort);

        System.out.println(" dt " + dt);
        System.out.println(" gV " + params.gV);
        System.out.println(" gVV " + params.gVV);

        //Ulat
        sim.integrator.reset();//uLat
        MeterPotentialEnergyFromIntegrator meterPE = new MeterPotentialEnergyFromIntegrator(sim.integrator);
        double ULat = meterPE.getDataAsScalar();
        System.out.println(" uLat " + (ULat/N));

        //Plat
        MeterPressure meterP = new MeterPressure(sim.box, sim.potentialMaster);
        meterP.setTemperature(0);//ZERO means NO ideal gas component!
        double pLat = meterP.getDataAsScalar();
        System.out.println(" pLat " + pLat);
        meterP.setTemperature(T);

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



        MeterSolidHMAFD meterElastic = new MeterSolidHMAFD(sim.space, new MeterPotentialEnergyFromIntegrator(sim.integrator), sim.potentialMaster,sim.potentialMasterp,sim.potentialMasterm, sim.boxp, sim.boxm, sim.coordinateDefinition, T, elasticParams, dbeta, deV);

//Equilibration
        final long startTime = System.currentTimeMillis();
        System.out.flush();
        sim.initialize(stepsEq);

        int numBlocks = 100;
        int interval = 10;
        long blockSize = steps / (numBlocks * interval);
        if (blockSize == 0) blockSize = 1;
        System.out.println(" block size " + blockSize + " interval " + interval);


        //Mapped
        AccumulatorAverageCovariance accumulatorElastic = new AccumulatorAverageCovariance(blockSize);
        DataPumpListener accumulatorUPPump = new DataPumpListener(meterElastic, accumulatorElastic, interval);
        sim.integrator.getEventManager().addListener(accumulatorUPPump);


        // Short run for P1_avg: to solve the Var(x) issue
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, stepsShort));
        IData dataElasticAvg0 = accumulatorElastic.getData(accumulatorElastic.AVERAGE);
        double uShift = dataElasticAvg0.getValue(1); //U_hma
        double pShift = dataElasticAvg0.getValue(3); //P_hma
        System.out.println(" uShift " + uShift);
        System.out.println(" pShift " + pShift);
        meterElastic.setShift(uShift, pShift);
        accumulatorElastic.reset();


// Production run ...
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, steps));

        IData dataElasticAvg = accumulatorElastic.getData(accumulatorElastic.AVERAGE);
        IData dataElasticErr = accumulatorElastic.getData(accumulatorElastic.ERROR);
        IData dataElasticCorr = accumulatorElastic.getData(accumulatorElastic.BLOCK_CORRELATION);
        IData dataElasticCov = accumulatorElastic.getData(accumulatorElastic.COVARIANCE);


/** First Derivatives*/
        double U_conv     = dataElasticAvg.getValue(0) + uShift;
        double errU_conv  = dataElasticErr.getValue(0);
        double corrU_conv = dataElasticCorr.getValue(0);

        double U_hma      = dataElasticAvg.getValue(1) + uShift;
        double errU_hma   = dataElasticErr.getValue(1);
        double corrU_hma  = dataElasticCorr.getValue(1);

        double P_conv     = dataElasticAvg.getValue(2) + pShift;
        double errP_conv  = dataElasticErr.getValue(2);
        double corrP_conv = dataElasticCorr.getValue(2);

        double P_hma      = dataElasticAvg.getValue(3) + pShift;
        double errP_hma   = dataElasticErr.getValue(3);
        double corrP_hma  = dataElasticCorr.getValue(3);


        System.out.println("\n++++ First Derivatives ++++");
        System.out.println(" U_conv   " + U_conv + " +/- " + errU_conv + "   corr " + corrU_conv);
        System.out.println(" U_hma    " + U_hma + " +/- " + errU_hma + "   corr " + corrU_hma);
        System.out.println(" P_conv   " + P_conv + " +/- " + errP_conv + "   corr " + corrP_conv);
        System.out.println(" P_hma    " + P_hma + " +/- " + errP_hma + "   corr " + corrP_hma);


        /** Second Derivatives*/
        int nd = dataElasticAvg.getLength();
        double varU_conv  = dataElasticCov.getValue(0 * (nd + 1));
        double varU_hma   = dataElasticCov.getValue(1 * (nd + 1));
        double varP_conv  = dataElasticCov.getValue(2 * (nd + 1));
        double varP_hma   = dataElasticCov.getValue(3 * (nd + 1));
        double T2 = T * T;
        // Cv = - d2A/db2 / T2 = [d2J - d2U + var] / T2
        double Cv_conv = N * varU_conv / T2;
        double Cv_hma  = N * varU_hma / T2 - dataElasticAvg.getValue(4) / T2 + 3.0 / 2.0 * (N - 1.0) / N ;

        double B_conv  =            rho * T   + dataElasticAvg.getValue(5) - V / T * varP_conv;
        double B_hma   = -3.0 * rho * T * gVV + dataElasticAvg.getValue(6) - V / T * varP_hma;

        System.out.println("\n++++ Second Derivatives ++++");
        System.out.println(" Cv_conv  " + Cv_conv);
        System.out.println(" Cv_hma   " + Cv_hma);

        System.out.println(" B_conv   " + B_conv);
        System.out.println(" B_hma    " + B_hma);


        long endTime = System.currentTimeMillis();
        System.out.println("\n Time " + (endTime - startTime) / 1000.0);
    }

    public static class SimLJPropsCijParam extends ParameterBase {
        public int N = 256;
        public double rho = 1;
        public long steps = 10000;
        public long stepsEq = steps / 10;
        public double T = 0.5;
        public double rc = 3;
        public double dt = 0.005;
        public double gV = 3.096385101500418;
        public double gVV = -0.200069081433023;
        public double dbeta = 0.001;
        public double deV = 0.001;
        public boolean isGraphic = !true;
        public boolean isLRC = false;
    }
}