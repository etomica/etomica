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
import etomica.data.types.DataGroup;
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

    public SimLJPropsCijFD(Space _space, int nMol, double density, double temperature, double rc, double dt, double deV) {
        super(_space);
        species = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this), true);
        addSpecies(species);

        double L = Math.pow(4.0/density, 1.0 / 3.0);
        int n = (int) Math.round(Math.pow(nMol / 4, 1.0 / 3.0));
        boundary = new BoundaryRectangularPeriodic(space, n * L);
        box = this.makeBox(boundary);
        box.setNMolecules(species, nMol);
        primitive = new PrimitiveCubic(space, n * L);
        nCells = new int[]{n, n, n};
        Basis basisFCC = new BasisCubicFcc();
        basis = new BasisBigCell(space, basisFCC, nCells);
        coordinateDefinition = new CoordinateDefinitionLeaf(box, primitive, basis, space);
        coordinateDefinition.initializeCoordinates(new int[]{1, 1, 1});
        potentialMaster = new PotentialMasterList(getSpeciesManager(), box, 2, rc, BondingInfo.noBonding());
        potentialMaster.doAllTruncationCorrection = false;
        integrator = new IntegratorVelocityVerlet(potentialMaster, getRandom(), dt, temperature, box);
        integrator.setIsothermal(true);
        integrator.setThermostatNoDrift(true);
        IPotential2 potential = new P2LennardJones();
        potential = new P2SoftSphericalTruncated(potential, rc);
        AtomType sphereType = species.getLeafType();
        potentialMaster.setPairPotential(sphereType, sphereType, potential);
        potentialMaster.init(); // should not need to do this in MeterSolidHMAFD
        integrator.getEventManager().removeListener(potentialMaster);
        this.getController().addActivity(new ActivityIntegrate(integrator));
        ((P2SoftSphericalTruncated) potential).setTruncationRadius(0.6 * boundary.getBoxSize().getX(0));

        //box+
        double densityp = density / (1 + deV);
        double Lp = Math.pow(4.0/densityp, 1.0 / 3.0);
        boundaryp = new BoundaryRectangularPeriodic(space, n * Lp);
        boxp = this.makeBox(boundaryp);
        boxp.setNMolecules(species, nMol);
        primitivep = new PrimitiveCubic(space, n * Lp);
        Basis basisFCCp = new BasisCubicFcc();
        basisp = new BasisBigCell(space, basisFCCp, nCells);
        coordinateDefinitionp = new CoordinateDefinitionLeaf(boxp, primitivep, basisp, space);
        coordinateDefinitionp.initializeCoordinates(new int[]{1, 1, 1});
        potentialMasterp = new PotentialMasterList(getSpeciesManager(), boxp, 2, rc, BondingInfo.noBonding());
        potentialMasterp.doAllTruncationCorrection = false;
        potentialMasterp.setPairPotential(sphereType, sphereType, potential);
        potentialMasterp.init();

        //box-
        double densitym = density / (1 - deV);
        double Lm = Math.pow(4.0/densitym, 1.0 / 3.0);
        boundarym = new BoundaryRectangularPeriodic(space, n * Lm);
        boxm = this.makeBox(boundarym);
        boxm.setNMolecules(species, nMol);
        primitivem = new PrimitiveCubic(space, n * Lm);
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
        long numSteps = params.numSteps;
        int nMol = params.nMol;
        double density = params.density;
        double rc = params.rc;
        double dbeta = params.dbeta;
        double deV = params.deV;
        double temperature = params.temperature;
        double dt = params.dt;
        boolean isGraphic = params.isGraphic;
        double[] elasticParams = new double[11];
        elasticParams[0] = params.gV;
        elasticParams[1] = params.gVV;
        System.out.println(" gV: " + params.gV + " gVV: " + params.gVV);
        System.out.println();


        final SimLJPropsCijFD sim = new SimLJPropsCijFD(Space.getInstance(3), nMol, density, temperature, rc, dt, deV);
        sim.integrator.reset();//uLat

        System.out.println(" " + nMol + " atoms " + " temperature " + temperature + " rc " + rc);
        System.out.println(" " + numSteps + " MD steps");
        System.out.println(" dt " + dt);
        double volume = sim.box.getBoundary().volume();
        System.out.println(" density " + density + " density " + density + " volume " + volume + "\n");

        MeterPotentialEnergyFromIntegrator meterPE = new MeterPotentialEnergyFromIntegrator(sim.integrator);
        double ULat = meterPE.getDataAsScalar();
        System.out.println(" uLat: " + (ULat/nMol));
        //P
        MeterPressure meterP = new MeterPressure(sim.box, sim.potentialMaster);
        meterP.setTemperature(0);//ZERO means NO ideal gas component!
        double pLat = meterP.getDataAsScalar();
        meterP.setTemperature(temperature);
        System.out.println(" pLat: " + pLat);


        MeterSolidHMAFD meterElastic = new MeterSolidHMAFD(sim.space, new MeterPotentialEnergyFromIntegrator(sim.integrator), sim.potentialMaster,sim.potentialMasterp,sim.potentialMasterm, sim.boxp, sim.boxm, sim.coordinateDefinition, temperature, elasticParams, dbeta, deV);

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

            //add
            simGraphic.add(historyPE);
            simGraphic.add(historyP);
            return;
        }

//Equilibration
        final long startTime = System.currentTimeMillis();
        System.out.flush();
        sim.initialize(numSteps / 10);

        int numBlocks = 100;
        int interval = 10;
        long blockSize = numSteps / (numBlocks * interval);
        if (blockSize == 0) blockSize = 1;
        System.out.println(" block size " + blockSize + " interval " + interval);

        //U
//        AccumulatorAverageFixed accumulatorPE = new AccumulatorAverageFixed(blockSize);
//        DataPumpListener accumulatorPEPump = new DataPumpListener(meterPE, accumulatorPE, interval);
//        sim.integrator.getEventManager().addListener(accumulatorPEPump);
//
//        //P
//        AccumulatorAverageFixed accumulatorP = new AccumulatorAverageFixed(blockSize);
//        DataPumpListener accumulatorPPump = new DataPumpListener(meterP, accumulatorP, interval);
//        sim.integrator.getEventManager().addListener(accumulatorPPump);


        //Mapped
        AccumulatorAverageCovariance accumulatorElastic = new AccumulatorAverageCovariance(blockSize);
        DataPumpListener accumulatorUPPump = new DataPumpListener(meterElastic, accumulatorElastic, interval);
        sim.integrator.getEventManager().addListener(accumulatorUPPump);


// Production run ...
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, numSteps));

        DataGroup dataElastic = (DataGroup) accumulatorElastic.getData();
        IData dataElasticAvg  = dataElastic.getData(accumulatorElastic.AVERAGE.index);
        IData dataElasticErr  = dataElastic.getData(accumulatorElastic.ERROR.index);
        IData dataElasticCorr = dataElastic.getData(accumulatorElastic.BLOCK_CORRELATION.index);

/** First Derivatives*/
        //Bulk
        double U_conv     = dataElasticAvg.getValue(0);
        double errU_conv  = dataElasticErr.getValue(0);
        double corrU_conv = dataElasticCorr.getValue(0);
        double U_hma      = dataElasticAvg.getValue(1);
        double errU_hma   = dataElasticErr.getValue(1);
        double corrU_hma  = dataElasticCorr.getValue(1);
        double P_conv     = dataElasticAvg.getValue(2);
        double errP_conv  = dataElasticErr.getValue(2);
        double corrP_conv = dataElasticCorr.getValue(2);
        double P_hma      = dataElasticAvg.getValue(3);
        double errP_hma   = dataElasticErr.getValue(3);
        double corrP_hma  = dataElasticCorr.getValue(3);

        System.out.println("\n++++ First Derivatives ++++");
        System.out.println(" U_conv   " + U_conv + " +/- " + errU_conv + "   corr: " + corrU_conv);
        System.out.println(" U_hma    " + U_hma + " +/- " + errU_hma + "   corr: " + corrU_hma);
        System.out.println(" P_conv   " + P_conv + " +/- " + errP_conv + "   corr: " + corrP_conv);
        System.out.println(" P_hma    " + P_hma + " +/- " + errP_hma + "   corr: " + corrP_hma);

        System.out.println();

        long endTime = System.currentTimeMillis();
        System.out.println("\n Time: " + (endTime - startTime) / 1000.0);
    }

    public static class SimLJPropsCijParam extends ParameterBase {
        public int nMol = 256;
        public double density = 1;
        public long numSteps = 100000;
        public double temperature = 0.1;
        public boolean isGraphic = !true;
        public double rc = 3;
        public double dt = 0.005;
        public double gV = 3.096385101500418;
        public double gVV = -0.200069081433023;
        public double dbeta = 0.001;
        public double deV = 0.00001;

    }
}