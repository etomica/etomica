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
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.lattice.crystal.*;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.*;
import etomica.simulation.Simulation;
import etomica.space.*;
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
    public int nMol;
    public Primitive primitive;
    public MCMoveAtomCoupled atomMove;
    public PotentialMasterList potentialMaster;
    public Potential2SoftSpherical potential;
    public SpeciesGeneral species;


    public SimLJPropsCij(Space _space, int nMol, double density0, double temperature, double rc, boolean isLRC, double strain_x, double strain_yz, boolean isSS, boolean isBCC, int nExp) {
        super(_space);
//        setRandom(new RandomMersenneTwister(1)); // set seed
        species = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this));
        addSpecies(species);


        int nBasis = isBCC ? 2 : 4;
        int nc = (int)Math.round(Math.pow(nMol/nBasis, 1.0/3.0));
        double L = Math.pow(nBasis / density0, 1.0 / 3.0);

        System.out.println(" rc = " + rc);

        nCells = new int[]{nc, nc, nc};
        if(isBCC){
            basis = new BasisCubicBcc();
        }else{
            basis = new BasisCubicFcc();
        }
        Vector[] cellDim = new Vector[3];
        cellDim[0] = Vector.of(new double[]{L, 0, 0});
        cellDim[1] = Vector.of(new double[]{0, L, 0});
        cellDim[2] = Vector.of(new double[]{0, 0, L});
        primitive = new PrimitiveGeneral(space, cellDim);
        boundary = new BoundaryDeformableLattice(primitive, nCells);
        box = this.makeBox(boundary);
        box.setNMolecules(species, nMol);

        CoordinateDefinitionLeaf coordinateDefinition0 = new CoordinateDefinitionLeaf(box, primitive, basis, space);
        coordinateDefinition0.initializeCoordinates(nCells);

        potentialMaster = new PotentialMasterList(this, space);
        integrator = new IntegratorMC(potentialMaster, getRandom(), temperature, box);
        MeterPotentialEnergy meterPE = new MeterPotentialEnergy(potentialMaster, box);
        atomMove = new MCMoveAtomCoupled(potentialMaster, meterPE, getRandom(), space);
        atomMove.setStepSize(0.1);
        atomMove.setStepSizeMax(0.5);
        atomMove.setDoExcludeNonNeighbors(true);
        integrator.getMoveManager().addMCMove(atomMove);

        Potential2SoftSpherical potential;
        if (isSS) {
            potential = new P2SoftSphere(space, 1.0, 1.0, nExp);
            System.out.println("** SS potential **");
        } else {
            potential = new P2LennardJones(space);
            System.out.println("** LJ potential **");
        }
        potential = new P2SoftSphericalTruncated(space, potential, rc);

        atomMove.setPotential(potential);

        AtomType sphereType = species.getLeafType();
        potentialMaster.addPotential(potential, new AtomType[]{sphereType, sphereType});

        potentialMaster.lrcMaster().setEnabled(isLRC);

        int cellRange = 4;
        potentialMaster.setRange(rc);
        potentialMaster.setCellRange(cellRange);
        potentialMaster.getNeighborManager(box).reset();
        int potentialCells = potentialMaster.getNbrCellManager(box).getLattice().getSize()[0];
        if (potentialCells < cellRange * 2 + 1) {
            throw new RuntimeException("oops (" + potentialCells + " < " + (cellRange * 2 + 1) + ")");
        }

        getController().addActivity(new ActivityIntegrate(integrator));

        if (!isLRC) {
            ((P2SoftSphericalTruncated) potential).setTruncationRadius(0.6 * boundary.getBoxSize().getX(0));
        }

        if (strain_x != 0 || strain_yz != 0) {
            cellDim[0] = Vector.of(new double[]{(1 + strain_x) * nc * L, 0, 0});
            cellDim[1] = Vector.of(new double[]{0, nc * L, strain_yz / 2 * nc * L});
            cellDim[2] = Vector.of(new double[]{0, strain_yz / 2 * nc * L, nc * L});

            ((BoundaryDeformablePeriodic) boundary).setEdgeVector(0, cellDim[0]);
            ((BoundaryDeformablePeriodic) boundary).setEdgeVector(1, cellDim[1]);
            ((BoundaryDeformablePeriodic) boundary).setEdgeVector(2, cellDim[2]);

            cellDim[0] = Vector.of(new double[]{(1 + strain_x) * L, 0, 0});
            cellDim[1] = Vector.of(new double[]{0, L, strain_yz / 2 * L});
            cellDim[2] = Vector.of(new double[]{0, strain_yz / 2 * L, L});

            primitive = new PrimitiveGeneral(space, cellDim);
            coordinateDefinition = new CoordinateDefinitionLeaf(box, primitive, basis, space);
            coordinateDefinition.initializeCoordinates(new int[]{nc, nc, nc});

        } else {
            coordinateDefinition = coordinateDefinition0;
        }

    }
    // End of Constructor

    public void initialize(long initSteps) {
        this.getController().runActivityBlocking(new ActivityIntegrate(this.integrator, initSteps));
        integrator.getMoveManager().setEquilibrating(false);
    }


    public static void main(String[] args) {
        SimLJPropsCijParam params = new SimLJPropsCijParam();
        ParseArgs.doParseArgs(params, args);
        long numSteps = params.numSteps;
        int nMol = params.nMol;
        int nExp = params.nExp;
        double density0 = params.density0;
        double rc = params.rc;
        double temperature = params.temperature;
        boolean isLRC = params.isLRC;
        boolean isSS = params.isSS;
        boolean isBCC = params.isBCC;
        boolean isGraphic = params.isGraphic;
        double strain_x = params.strain_x;
        double strain_yz = params.strain_yz;
        double[] elasticParams = new double[11];
        elasticParams[0] = params.gV;
        elasticParams[1] = params.gVV;
        elasticParams[2] = params.gx1;
        elasticParams[3] = params.gy1;
        elasticParams[4] = params.gy4;
        elasticParams[5] = params.gx11;
        elasticParams[6] = params.gy11;
        elasticParams[7] = params.gx44;
        elasticParams[8] = params.gy44;
        elasticParams[9] = params.gx12;
        elasticParams[10] = params.gz12;
        System.out.println(" HMA Elastic Parameters:" + (isSS ? " SS-"+nExp : " LJ") + (isBCC ? " BCC " : " FCC"));
        System.out.println(" gV: " + params.gV + " gVV: " + params.gVV);
        System.out.println(" gx1: " + params.gx1 + " gy1: " + params.gy1 + " gy4 " + params.gy4);
        System.out.println(" gx11: " + params.gx11 + " gy11: " + params.gy11 + " gx44 " + params.gx44);
        System.out.println(" gy44: " + params.gy44 + " gx12: " + params.gx12 + " gz12 " + params.gz12);
        System.out.println();


        final SimLJPropsCij sim = new SimLJPropsCij(Space.getInstance(3), nMol, density0, temperature, rc, isLRC, strain_x, strain_yz, isSS, isBCC, nExp);
        sim.integrator.reset(); // so we can ask it for the PE (the integrator already gets reset at the start of the simulation)

        System.out.println(" " + nMol + " atoms " + " temperature " + temperature + " rc " + rc);
        System.out.println(" " + numSteps + " MC steps");
        System.out.println(" strain_x  " + strain_x + " strain_yz  " + strain_yz);

        double volume = sim.box.getBoundary().volume();
        double density = nMol/volume;
        System.out.println(" density0 " + density0 + " density " + density + " volume " + volume + "\n");

        MeterPotentialEnergyFromIntegrator meterPE = new MeterPotentialEnergyFromIntegrator(sim.integrator);
        double ULat = meterPE.getDataAsScalar();
        System.out.println(" u_lat " + (ULat/nMol));
        //P
        MeterPressure meterP = new MeterPressure(sim.space);
        meterP.setBox(sim.box);
        meterP.setTemperature(0);//ZERO means NO ideal gas component!
        meterP.setPotentialMaster(sim.potentialMaster);
        double PLat = meterP.getDataAsScalar();
        meterP.setTemperature(temperature);
        System.out.println(" P_lat " + PLat);

        //Pij
        MeterPressureTensor meterPij = new MeterPressureTensor(sim.potentialMaster, sim.space);
        meterPij.setBox(sim.box);
        meterPij.setTemperature(0);
//        meterPij.getPotentialCalculation().setIntegrator(sim.integrator);
        IData dataPij_lat = meterPij.getData();
        meterPij.setTemperature(temperature);
        double P1_lat = -dataPij_lat.getValue(0);
        double P4_lat = -dataPij_lat.getValue(5);
        System.out.println(" P1_lat " + P1_lat);
        System.out.println(" P4_lat " + P4_lat);

        //Elastic constants
        MeterSolidPropsLJ meterElasticLat = new MeterSolidPropsLJ(sim.space, new MeterPotentialEnergyFromIntegrator(sim.integrator), sim.potentialMaster, sim.coordinateDefinition, temperature, elasticParams);
        IData dataCijLat = meterElasticLat.getData();
        double B_lat = dataCijLat.getValue(17) - density * temperature;
        double bV_lat = dataCijLat.getValue(37) - density * temperature;
        double C11_lat = dataCijLat.getValue(19) - 2 * density * temperature;
        double C12_lat = dataCijLat.getValue(25);
        double C44_lat = dataCijLat.getValue(31) - density * temperature;

        System.out.println(" B_lat   " + B_lat);
        System.out.println(" C11_lat " + C11_lat);
        System.out.println(" C12_lat " + C12_lat);
        System.out.println(" C44_lat " + C44_lat);
        System.out.println();

        MeterSolidPropsLJ meterElastic = new MeterSolidPropsLJ(sim.space, new MeterPotentialEnergyFromIntegrator(sim.integrator), sim.potentialMaster, sim.coordinateDefinition, temperature, elasticParams);

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
            simGraphic.makeAndDisplayFrame(("LJ-") + (isBCC ? "BCC" : " FCC"));
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
            // P1
            DisplayPlot historyPij = new DisplayPlot();
            historyPij.setLabel("Pij");
            // Conv P1
            AccumulatorHistory accPij = new AccumulatorHistory(new HistoryCollapsingAverage());
            accPij.setTimeDataSource(timeSource);
            DataSplitter splitter = new DataSplitter();
            DataPumpListener pumpPij = new DataPumpListener(meterPij, splitter, 10);
            splitter.setDataSink(0, accPij);
            sim.integrator.getEventManager().addListener(pumpPij);
            accPij.setDataSink(historyPij.getDataSet().makeDataSink());
            historyPij.setLegend(new DataTag[]{accPij.getTag()}, "Conv");
            // Mapped P1
            DataSplitter splitterM = new DataSplitter();
            DataPumpListener pumpPijM = new DataPumpListener(meterElastic, splitterM, 10);
            //add
            simGraphic.add(historyPE);
            simGraphic.add(historyP);
            simGraphic.add(historyPij);
            return;
        }

//Equilibration
        final long startTime = System.currentTimeMillis();
        System.out.flush();
        sim.initialize(numSteps / 5); //20%


        int numBlocks = 100;
        int interval = nMol;
        long blockSize = numSteps / (numBlocks * interval);
        if (blockSize == 0) blockSize = 1;
        System.out.println(" block size " + blockSize + " interval " + interval);
        //U
        AccumulatorAverageFixed accumulatorPE = new AccumulatorAverageFixed(blockSize);
        DataPumpListener accumulatorPEPump = new DataPumpListener(meterPE, accumulatorPE, interval);
        sim.integrator.getEventManager().addListener(accumulatorPEPump);
        //P
        AccumulatorAverageFixed accumulatorP = new AccumulatorAverageFixed(blockSize);
        DataPumpListener accumulatorPPump = new DataPumpListener(meterP, accumulatorP, interval);
        sim.integrator.getEventManager().addListener(accumulatorPPump);

        //Pij
        AccumulatorAverageFixed accumulatorPij = new AccumulatorAverageFixed(blockSize);
        DataPumpListener accumulatorPijPump = new DataPumpListener(meterPij, accumulatorPij, interval);
        sim.integrator.getEventManager().addListener(accumulatorPijPump);

        //Mapped
        AccumulatorAverageCovariance accumulatorElastic = new AccumulatorAverageCovariance(blockSize);
        DataPumpListener accumulatorUPPump = new DataPumpListener(meterElastic, accumulatorElastic, interval);
        sim.integrator.getEventManager().addListener(accumulatorUPPump);

// Short run for P1_avg: to solve the Var(x) issue
        long Nshort = numSteps / 10;
        System.out.println(" N_short_sim = " + Nshort);
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, Nshort));
        DataGroup dataElastic0 = (DataGroup) accumulatorElastic.getData();
        IData dataElasticAvg0 = dataElastic0.getData(accumulatorElastic.AVERAGE.index);
        double Ushift = dataElasticAvg0.getValue(1); //U_hma
        double Pshift = dataElasticAvg0.getValue(3); //P_hma
        System.out.println(" Ushift: " + Ushift + "   Pshift: " + Pshift);
        meterElastic.setShift(Ushift, Pshift);
        accumulatorElastic.reset();


// Production run ...
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, numSteps));

        DataGroup dataElastic = (DataGroup) accumulatorElastic.getData();
        IData dataElasticAvg  = dataElastic.getData(accumulatorElastic.AVERAGE.index);
        IData dataElasticErr  = dataElastic.getData(accumulatorElastic.ERROR.index);
        IData dataElasticCorr = dataElastic.getData(accumulatorElastic.BLOCK_CORRELATION.index);
        IData dataElasticCov  = dataElastic.getData(accumulatorElastic.COVARIANCE.index);

/** First Derivatives*/
        //Bulk
        double U_conv     = dataElasticAvg.getValue(0) + Ushift;
        double errU_conv  = dataElasticErr.getValue(0);
        double corrU_conv = dataElasticCorr.getValue(0);
        double U_hma      = dataElasticAvg.getValue(1) + Ushift;
        double errU_hma   = dataElasticErr.getValue(1);
        double corrU_hma  = dataElasticCorr.getValue(1);
        double P_conv     = dataElasticAvg.getValue(2) + Pshift;
        double errP_conv  = dataElasticErr.getValue(2);
        double corrP_conv = dataElasticCorr.getValue(2);
        double P_hma      = dataElasticAvg.getValue(3) + Pshift;
        double errP_hma   = dataElasticErr.getValue(3);
        double corrP_hma  = dataElasticCorr.getValue(3);
        //Normal stress
        double P1_conv     = dataElasticAvg.getValue(4) - Pshift;
        double errP1_conv  = dataElasticErr.getValue(4);
        double corrP1_conv = dataElasticCorr.getValue(4);
        double P1_hma      = dataElasticAvg.getValue(5) - Pshift;
        double errP1_hma   = dataElasticErr.getValue(5);
        double corrP1_hma  = dataElasticCorr.getValue(5);
        double P2_conv     = dataElasticAvg.getValue(6) - Pshift;
        double errP2_conv  = dataElasticErr.getValue(6);
        double corrP2_conv = dataElasticCorr.getValue(6);
        double P2_hma      = dataElasticAvg.getValue(7) - Pshift;
        double errP2_hma   = dataElasticErr.getValue(7);
        double corrP2_hma  = dataElasticCorr.getValue(7);
        double P3_conv     = dataElasticAvg.getValue(8) - Pshift;
        double errP3_conv  = dataElasticErr.getValue(8);
        double corrP3_conv = dataElasticCorr.getValue(8);
        double P3_hma      = dataElasticAvg.getValue(9) - Pshift;
        double errP3_hma   = dataElasticErr.getValue(9);
        double corrP3_hma  = dataElasticCorr.getValue(9);
        //Shear stress
        double P4_conv     = dataElasticAvg.getValue(10);
        double errP4_conv  = dataElasticErr.getValue(10);
        double corrP4_conv = dataElasticCorr.getValue(10);
        double P4_hma      = dataElasticAvg.getValue(11);
        double errP4_hma   = dataElasticErr.getValue(11);
        double corrP4_hma  = dataElasticCorr.getValue(11);
        double P5_conv     = dataElasticAvg.getValue(12);
        double errP5_conv  = dataElasticErr.getValue(12);
        double corrP5_conv = dataElasticCorr.getValue(12);
        double P5_hma      = dataElasticAvg.getValue(13);
        double errP5_hma   = dataElasticErr.getValue(13);
        double corrP5_hma  = dataElasticCorr.getValue(13);
        double P6_conv     = dataElasticAvg.getValue(14);
        double errP6_conv  = dataElasticErr.getValue(14);
        double corrP6_conv = dataElasticCorr.getValue(14);
        double P6_hma      = dataElasticAvg.getValue(15);
        double errP6_hma   = dataElasticErr.getValue(15);
        double corrP6_hma  = dataElasticCorr.getValue(15);

/** Second Derivatives*/
        int nd = dataElasticAvg.getLength();
        double varU_conv  = dataElasticCov.getValue(0 * (nd + 1));
        double varU_hma   = dataElasticCov.getValue(1 * (nd + 1));
        double varP_conv  = dataElasticCov.getValue(2 * (nd + 1));
        double varP_hma   = dataElasticCov.getValue(3 * (nd + 1));
        double varP1_conv = dataElasticCov.getValue(4 * (nd + 1));
        double varP1_hma  = dataElasticCov.getValue(5 * (nd + 1));
        double varP2_conv = dataElasticCov.getValue(6 * (nd + 1));
        double varP2_hma  = dataElasticCov.getValue(7 * (nd + 1));
        double varP3_conv = dataElasticCov.getValue(8 * (nd + 1));
        double varP3_hma  = dataElasticCov.getValue(9 * (nd + 1));
        double varP4_conv = dataElasticCov.getValue(10 * (nd + 1));
        double varP4_hma  = dataElasticCov.getValue(11 * (nd + 1));
        double varP5_conv = dataElasticCov.getValue(12 * (nd + 1));
        double varP5_hma  = dataElasticCov.getValue(13 * (nd + 1));
        double varP6_conv = dataElasticCov.getValue(14 * (nd + 1));
        double varP6_hma  = dataElasticCov.getValue(15 * (nd + 1));

        double covP1P2_conv = dataElasticCov.getValue(4 * nd + 6);
        double covP1P2_hma  = dataElasticCov.getValue(5 * nd + 7);
        double covP1P3_conv = dataElasticCov.getValue(4 * nd + 8);
        double covP1P3_hma  = dataElasticCov.getValue(5 * nd + 9);
        double covP2P3_conv = dataElasticCov.getValue(6 * nd + 8);
        double covP2P3_hma  = dataElasticCov.getValue(7 * nd + 9);

        double covUP_conv  = dataElasticCov.getValue(0 * nd + 2);
        double covUP_hma   = dataElasticCov.getValue(1 * nd + 3);
        double covUP1_conv = dataElasticCov.getValue(0 * nd + 4);
        double covUP1_hma  = dataElasticCov.getValue(1 * nd + 5);
        double covUP2_conv = dataElasticCov.getValue(0 * nd + 6);
        double covUP2_hma  = dataElasticCov.getValue(1 * nd + 7);
        double covUP3_conv = dataElasticCov.getValue(0 * nd + 8);
        double covUP3_hma  = dataElasticCov.getValue(1 * nd + 9);


//Bulk
        double Cv_conv = varU_conv / temperature / temperature + 3.0/2.0*(nMol-1);
        double Cv_hma  = dataElasticAvg.getValue(16) + varU_hma / temperature / temperature + 3.0*(nMol-1.0);
        double B_conv  = dataElasticAvg.getValue(17) - volume / temperature * varP_conv;
        double B_hma   = dataElasticAvg.getValue(18) - volume / temperature * varP_hma;
//Cij
        //normal
        double C11_conv = dataElasticAvg.getValue(19) - volume / temperature * varP1_conv;
        double C11_hma  = dataElasticAvg.getValue(20) - volume / temperature * varP1_hma;
        double C22_conv = dataElasticAvg.getValue(21) - volume / temperature * varP2_conv;
        double C22_hma  = dataElasticAvg.getValue(22) - volume / temperature * varP2_hma;
        double C33_conv = dataElasticAvg.getValue(23) - volume / temperature * varP3_conv;
        double C33_hma  = dataElasticAvg.getValue(24) - volume / temperature * varP3_hma;

        double C12_conv = dataElasticAvg.getValue(25) - volume / temperature * covP1P2_conv;
        double C12_hma  = dataElasticAvg.getValue(26) - volume / temperature * covP1P2_hma;
        double C13_conv = dataElasticAvg.getValue(27) - volume / temperature * covP1P3_conv;
        double C13_hma  = dataElasticAvg.getValue(28) - volume / temperature * covP1P3_hma;
        double C23_conv = dataElasticAvg.getValue(29) - volume / temperature * covP2P3_conv;
        double C23_hma  = dataElasticAvg.getValue(30) - volume / temperature * covP2P3_hma;
        //shear
        double C44_conv = dataElasticAvg.getValue(31) - volume / temperature * varP4_conv;
        double C44_hma  = dataElasticAvg.getValue(32) - volume / temperature * varP4_hma;
        double C55_conv = dataElasticAvg.getValue(33) - volume / temperature * varP5_conv;
        double C55_hma  = dataElasticAvg.getValue(34) - volume / temperature * varP5_hma;
        double C66_conv = dataElasticAvg.getValue(35) - volume / temperature * varP6_conv;
        double C66_hma  = dataElasticAvg.getValue(36) - volume / temperature * varP6_hma;

//b_mn
        double gV_conv  = density + 1.0 / temperature / temperature * covUP_conv;
        double gV_hma   = dataElasticAvg.getValue(37) + 1.0 / temperature / temperature * covUP_hma;
        double b11_conv = -density + 1.0 / temperature / temperature * covUP1_conv;
        double b11_hma  = dataElasticAvg.getValue(38) + 1.0 / temperature / temperature * covUP1_hma;
        double b22_conv = -density + 1.0 / temperature / temperature * covUP2_conv;
        double b22_hma  = dataElasticAvg.getValue(39) + 1.0 / temperature / temperature * covUP2_hma;
        double b33_conv = -density + 1.0 / temperature / temperature * covUP3_conv;
        double b33_hma  = dataElasticAvg.getValue(40) + 1.0 / temperature / temperature * covUP3_hma;


        System.out.println("\n++++ First Derivatives ++++");
        System.out.println(" U_conv   " + U_conv + " +/- " + errU_conv + "   corr: " + corrU_conv);
        System.out.println(" U_hma    " + U_hma + " +/- " + errU_hma + "   corr: " + corrU_hma);
        System.out.println(" P_conv   " + P_conv + " +/- " + errP_conv + "   corr: " + corrP_conv);
        System.out.println(" P_hma    " + P_hma + " +/- " + errP_hma + "   corr: " + corrP_hma);
        System.out.println(" P1_conv  " + P1_conv + " +/- " + errP1_conv + "   corr: " + corrP1_conv);
        System.out.println(" P1_hma   " + P1_hma + " +/- " + errP1_hma + "   corr: " + corrP1_hma);
        System.out.println(" P2_conv  " + P2_conv + " +/- " + errP2_conv + "   corr: " + corrP2_conv);
        System.out.println(" P2_hma   " + P2_hma + " +/- " + errP2_hma + "   corr: " + corrP2_hma);
        System.out.println(" P4_conv  " + P4_conv + " +/- " + errP4_conv + "   corr: " + corrP4_conv);
        System.out.println(" P4_hma   " + P4_hma + " +/- " + errP4_hma + "   corr: " + corrP4_hma);

        System.out.println("\n++++ Second Derivatives ++++");
        System.out.println(" Cv_conv  " + Cv_conv);
        System.out.println(" Cv_hma   " + Cv_hma);
        System.out.println(" B_conv   " + B_conv);
        System.out.println(" B_hma    " + B_hma);
        System.out.println(" gV_conv " + gV_conv);
        System.out.println(" gV_hma  " + gV_hma);
        System.out.println(" alphaV_conv " + gV_conv/B_conv);
        System.out.println(" alphaV_hma  " + gV_hma/B_hma);
        System.out.println(" Bs_conv  " + (B_conv+temperature*volume*gV_conv*gV_conv/Cv_conv));
        System.out.println(" Bs_hma   " + (B_hma+temperature*volume*gV_hma*gV_hma/Cv_hma));

        System.out.println();
        System.out.println(" C11_conv     " + C11_conv);
        System.out.println(" C11_conv_avg " + 1.0/3.0*(C11_conv + C22_conv + C33_conv));
        System.out.println(" C11_hma      " + C11_hma);
        System.out.println(" C11_hma_avg  " + 1.0/3.0*(C11_hma + C22_hma + C33_hma));
        System.out.println();
        System.out.println(" C12_conv     " + C12_conv);
        System.out.println(" C12_conv_avg " + 1.0/3.0*(C12_conv + C13_conv + C23_conv));
        System.out.println(" C12_hma      " + C12_hma);
        System.out.println(" C12_hma_avg  " + 1.0/3.0*(C12_hma + C13_hma + C23_hma));
        System.out.println();
        System.out.println(" C44_conv     " + C44_conv);
        System.out.println(" C44_conv_avg " + 1.0/3.0*(C44_conv + C55_conv + C66_conv));
        System.out.println(" C44_hma      " + C44_hma);
        System.out.println(" C44_hma_avg  " + 1.0/3.0*(C44_hma + C55_hma + C66_hma));
        System.out.println();
        System.out.println(" b11_conv     " + b11_conv);
        System.out.println(" b22_conv     " + b22_conv);
        System.out.println(" b11_conv_avg " + 1.0/3.0*(b11_conv + b22_conv + b33_conv));
        System.out.println(" b11_hma      " + b11_hma);
        System.out.println(" b22_hma      " + b22_hma);
        System.out.println(" b11_hma_avg  " + 1.0/3.0*(b11_hma + b22_hma + b33_hma));

    //Adiabatic Cij
        double C11s_conv = C11_conv + temperature*volume/Cv_conv*gV_conv*gV_conv;
        double C11s_hma  = C11_hma  + temperature*volume/Cv_hma*gV_hma*gV_hma;
        double C22s_conv = C22_conv + temperature*volume/Cv_conv*gV_conv*gV_conv;
        double C22s_hma  = C22_hma  + temperature*volume/Cv_hma*gV_hma*gV_hma;
        double C33s_conv = C33_conv + temperature*volume/Cv_conv*gV_conv*gV_conv;
        double C33s_hma  = C33_hma  + temperature*volume/Cv_hma*gV_hma*gV_hma;

        double C12s_conv = C12_conv + temperature*volume/Cv_conv*gV_conv*gV_conv;
        double C12s_hma  = C12_hma  + temperature*volume/Cv_hma*gV_hma*gV_hma;
        double C13s_conv = C13_conv + temperature*volume/Cv_conv*gV_conv*gV_conv;
        double C13s_hma  = C13_hma  + temperature*volume/Cv_hma*gV_hma*gV_hma;
        double C23s_conv = C23_conv + temperature*volume/Cv_conv*gV_conv*gV_conv;
        double C23s_hma  = C23_hma  + temperature*volume/Cv_hma*gV_hma*gV_hma;

        System.out.println("\n Adiabatic Cij");
        System.out.println(" C11s_conv     " + C11s_conv);
        System.out.println(" C11s_conv_avg " + 1.0/3.0*(C11s_conv + C22s_conv + C33s_conv));
        System.out.println(" C11s_hma      " + C11s_hma);
        System.out.println(" C11s_hma_avg  " + 1.0/3.0*(C11s_hma + C22s_hma + C33s_hma));
        System.out.println();
        System.out.println(" C12s_conv     " + C12s_conv);
        System.out.println(" C12s_conv_avg " + 1.0/3.0*(C12s_conv + C13s_conv + C23s_conv));
        System.out.println(" C12s_hma      " + C12_hma);
        System.out.println(" C12s_hma_avg  " + 1.0/3.0*(C12s_hma + C13s_hma + C23s_hma));
        System.out.println();

        long endTime = System.currentTimeMillis();
        System.out.println("\n Time: " + (endTime - startTime) / 1000.0);
    }

    public static class SimLJPropsCijParam extends ParameterBase {
        public boolean isBCC = true;
        public boolean isSS = true;
        public int nExp = 6;
        public int nMol = 500;
        public double density0 = 1.0;
        public long numSteps = 1000000;
        public double temperature = 0.1;
        public boolean isLRC = false;
        public boolean isGraphic = false;
        public double rc = 3.0;
        public double strain_x = 0.0;
        public double strain_yz = 0.0;

        //SS-FCC: nSS=12 , rho=1.0, rc=3.0, nC=5 (N=500)
//         public double gV=2.333333333332616 
//         public double gVV=-2.333333331008459
//         public double gx1=2.365809823734493
//         public double gy1=2.317095089949099
//         public double gy4=1.657754708346944
//         public double gx11=6.914084603034437
//         public double gy11=5.226926773132499
//         public double gx44=9.333235933502698
//         public double gy44=1.780559360999292
//         public double gx12=-5.822852205434982
//         public double gz12=-4.038264917932466


        //SS-BCC: nSS=6 , rho=1.0, rc=3.0, nC=5 (N=250)
        public double gV=1.333333333334733;
        public double gVV=-1.333333331020438;
        public double gx1=1.321732322767388;
        public double gy1=1.339133838679187;
        public double gy4=1.339133838700971;
        public double gx11=-5.208082731208532;
        public double gy11=-1.366754230260249;
        public double gx44=-6.180381361666975;
        public double gy44=1.328914959737209;
        public double gx12=1.282309044414773;
        public double gz12=-2.593822489037677;
    }
}
