package etomica.rotation;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomType;
import etomica.atom.DiameterHashByType;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.config.ConfigurationFile;
import etomica.config.ConfigurationFileBinary;
import etomica.data.*;
import etomica.data.history.HistoryCollapsingAverage;
import etomica.data.meter.MeterEnergy;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.exception.ConfigurationOverlapException;
import etomica.graphics.*;
import etomica.integrator.IntegratorMC;
import etomica.integrator.IntegratorVelocityVerlet.MyAgent;
import etomica.integrator.IntegratorVelocityVerletRattle;
import etomica.integrator.IntegratorVelocityVerletShake.BondConstraints;
import etomica.integrator.mcmove.MCMoveMolecule;
import etomica.integrator.mcmove.MCMoveRotateMolecule3D;
import etomica.listener.IntegratorListenerAction;
import etomica.models.clathrates.MinimizationTIP4P.ChargeAgentSourceRPM;
import etomica.models.water.ConformationWaterTIP4P;
import etomica.models.water.SpeciesWater4P;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.molecule.MoleculeAgentManager;
import etomica.molecule.MoleculePositionCOM;
import etomica.normalmode.*;
import etomica.potential.EwaldSummation;
import etomica.potential.P2LennardJones;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.IOrientation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.OrientationFull3D;
import etomica.space3d.Space3D;
import etomica.units.*;
import etomica.units.dimensions.Null;
import etomica.util.Constants;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.util.random.RandomMersenneTwister;

import java.awt.*;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.List;

/**
 * @author Weisong Lin
 */
public class WaterDADBMC extends Simulation {
    public IntegratorMC integrator;
    public PotentialMaster potentialMaster;
    public SpeciesWater4P species;
    public Box box;
    public ActivityIntegrate ai;
    protected Potential2SoftSphericalLS potentialLJLS;
    protected final MoleculeAgentManager latticeCoordinates;
    protected MCMoveRotateMolecule3D rotate;

    public WaterDADBMC(final Space space, double temperature, int numCells, double rCutRealES, double rCutLJ, boolean isIce, double kCut, boolean unitCells, final boolean doTranslation, final boolean doRotation) {
        super(space);
        setRandom(new RandomMersenneTwister(2));
        System.out.println("The random seed is set to be fixed");

        box = new Box(space);
        addBox(box);
        species = new SpeciesWater4P(getSpace(), false);
        addSpecies(species);
        box.setNMolecules(species, 46 * numCells * numCells * numCells);
        box.setDensity(46 / 12.03 / 12.03 / 12.03);
        ChargeAgentSourceRPM agentSource = new ChargeAgentSourceRPM(species, isIce);
        AtomLeafAgentManager<EwaldSummation.MyCharge> atomAgentManager = new AtomLeafAgentManager<EwaldSummation.MyCharge>(agentSource, box, EwaldSummation.MyCharge.class);
        int[] nC = new int[]{numCells, numCells, numCells};
        numCells = 1;
        ConfigurationFile config = new ConfigurationFile(numCells + "ncFinalPos");

        if (unitCells) {
            ConfigurationFileBinary.replicate(config, box, nC, space);
        } else {
            config.initializeCoordinates(box);
        }
        double a0 = box.getBoundary().getBoxSize().getX(0);
        double[] rC = new double[]{a0, a0, a0};
        double sigma, epsilon; //TIP4P
        if (isIce) {
            sigma = 3.1668;
            epsilon = Kelvin.UNIT.toSim(106.1);//TIP4P/Ice
        } else {//TIP4P
            double A = 600E3; // kcal A^12 / mol
            double C = 610.0; // kcal A^6 / mol
            double s6 = A / C;
            sigma = Math.pow(s6, 1.0 / 6.0);
            epsilon = Mole.UNIT.fromSim(Calorie.UNIT.toSim(C / s6 * 1000)) / 4.0;
        }
        latticeCoordinates = new MoleculeAgentManager(this, box, new MoleculeSiteSource(space, new MoleculePositionCOM(space), new WaterOrientationDefinition(space)));
        AtomLeafAgentManager<Vector> atomLatticeCoordinates = new AtomLeafAgentManager<>(new AtomSiteSource(space), box, Vector.class);
        EwaldSummationLattice potentialES = new EwaldSummationLattice(box, atomAgentManager, space, kCut, rCutRealES, atomLatticeCoordinates);
        P2LennardJones potentialLJ = new P2LennardJones(space, sigma, epsilon);
        potentialLJLS = new Potential2SoftSphericalLS(space, rCutLJ, rC, potentialLJ, latticeCoordinates);
//		potentialLJ =  new P2SoftSphericalTruncated(space, potentialLJ, rC);
        potentialMaster = new PotentialMaster();
        potentialMaster.addPotential(potentialES, new AtomType[0]);
        potentialMaster.addPotential(potentialLJLS, new AtomType[]{species.getOxygenType(), species.getOxygenType()});


        double lOH = ConformationWaterTIP4P.bondLengthOH;
        double lHH = Math.sqrt(2 * lOH * lOH * (1 - Math.cos(ConformationWaterTIP4P.angleHOH)));
        double lOM = ConformationWaterTIP4P.rOM;
        double lMH = Math.sqrt(lOH * lOH + lOM * lOM - 2 * lOH * lOM * Math.cos(0.5 * ConformationWaterTIP4P.angleHOH));

//        MCMoveMoleculeCoupled move = new MCMoveMoleculeCoupled(potentialMaster, getRandom(), space);
        MCMoveMolecule move = new MCMoveMolecule(this, potentialMaster, space);
        move.setBox(box);

        rotate = new MCMoveRotateMolecule3D(potentialMaster, getRandom(), space);
        rotate.setBox(box);

        integrator = new IntegratorMC(potentialMaster, getRandom(), Kelvin.UNIT.toSim(temperature));
        if (doTranslation) integrator.getMoveManager().addMCMove(move);
        if (doRotation) integrator.getMoveManager().addMCMove(rotate);

        integrator.setBox(box);
        integrator.setTemperature(Kelvin.UNIT.toSim(temperature));
        try {
            integrator.reset();
        } catch (ConfigurationOverlapException e) {
        }
        ai = new ActivityIntegrate(integrator);
        getController().addAction(ai);

    }

    public static void main(String[] args) {
        final long startTime = System.currentTimeMillis();
        WaterDADBParam waterDADBParam = new WaterDADBParam();
        ParseArgs.doParseArgs(waterDADBParam, args);
        final double temperature = waterDADBParam.temperature;
        int numCells = waterDADBParam.numCells;
        int numSteps = waterDADBParam.numSteps;
        double rCutRealES = waterDADBParam.rCutRealES;
        double rCutLJ = waterDADBParam.rCutLJ;
        boolean isIce = waterDADBParam.isIce;
        double kCut = waterDADBParam.kCut;
        boolean uniteCells = waterDADBParam.unitCells;
        boolean doTranslation = waterDADBParam.doTranslation;
        boolean doRotation = waterDADBParam.doRotation;
        boolean runGraphic = waterDADBParam.runGraphic;
        final WaterDADBMC sim = new WaterDADBMC(Space3D.getInstance(), temperature, numCells, rCutRealES, rCutLJ, isIce, kCut, uniteCells, doTranslation, doRotation);
        MeterPotentialEnergy meterPE2 = new MeterPotentialEnergy(sim.potentialMaster);
        meterPE2.setBox(sim.box);
        final double latticeEnergy = meterPE2.getDataAsScalar();
        System.out.println("latticeEnergy = " + latticeEnergy);

//      try{
//      	EwaldSummation.fileWriter.close();
//  	}
//  	catch (IOException e){
//  		throw new RuntimeException(e);
//  	}
//      System.exit(2);


        if (runGraphic) {
//			SpeciesSpheresMono guestSpecies = new SpeciesSpheresMono(sim,sim.space);
//			sim.addSpecies(guestSpecies);
//			sim.box.setNMolecules(guestSpecies, 8);
//        	sim.box.getMoleculeList(guestSpecies).getMolecule(0).getChildList().getAtom(0).getPosition().E(new double [] {6, 6, 6});
//        	sim.box.getMoleculeList(guestSpecies).getMolecule(1).getChildList().getAtom(0).getPosition().E(new double [] {6, -6, -6});
//        	sim.box.getMoleculeList(guestSpecies).getMolecule(2).getChildList().getAtom(0).getPosition().E(new double [] {6, -6, 6});
//        	sim.box.getMoleculeList(guestSpecies).getMolecule(3).getChildList().getAtom(0).getPosition().E(new double [] {6, 6, -6});
//        	sim.box.getMoleculeList(guestSpecies).getMolecule(4).getChildList().getAtom(0).getPosition().E(new double [] {-6, 6, 6});
//        	sim.box.getMoleculeList(guestSpecies).getMolecule(5).getChildList().getAtom(0).getPosition().E(new double [] {-6, -6, -6});
//        	sim.box.getMoleculeList(guestSpecies).getMolecule(6).getChildList().getAtom(0).getPosition().E(new double [] {-6, 6, -6});
//        	sim.box.getMoleculeList(guestSpecies).getMolecule(7).getChildList().getAtom(0).getPosition().E(new double [] {-6, -6, 6});
//			sim.box.getMoleculeList(guestSpecies).getMolecule(0).getChildList().getAtom(0).getPosition().E(new double [] {0, 0, 0});

//			sim.ai.setSleepPeriod(2);
            SimulationGraphic graphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, "Rattle", 1, sim.space, sim.getController());
            ((ColorSchemeByType) graphic.getDisplayBox(sim.box).getColorScheme()).setColor(sim.species.getHydrogenType(), Color.WHITE);
            ((ColorSchemeByType) graphic.getDisplayBox(sim.box).getColorScheme()).setColor(sim.species.getOxygenType(), Color.RED);
            ((DiameterHashByType) graphic.getDisplayBox(sim.box).getDiameterHash()).setDiameter(sim.species.getMType(), 0.1);
//			((DiameterHashByType)graphic.getDisplayBox(sim.box).getDiameterHash()).setDiameter(guestSpecies.getAtomType(0), 2.5);

//			((ColorSchemeByType)graphic.getDisplayBox(sim.box).getColorScheme()).setColor(guestSpecies.getAtomType(0), Color.ORANGE);

            ((DisplayBoxCanvasG3DSys) graphic.getDisplayBox(sim.box).canvas).setBackgroundColor(Color.WHITE);
            ((DisplayBoxCanvasG3DSys) graphic.getDisplayBox(sim.box).canvas).setBoundaryFrameColor(Color.black);
            //((DiameterHashByType)graphic.getDisplayBox(sim.box).getDiameterHash()).setDiameter(sim.species.getOxygenType(), .1);
//            ((DiameterHashByType)graphic.getDisplayBox(box).getDiameterHash()).setDiameter(hType, 1);
//            MeterEnergy meterE = new MeterEnergy(sim.potentialMaster, sim.box);
            AccumulatorHistory history = new AccumulatorHistory(new HistoryCollapsingAverage());
            history.setTimeDataSource(new DataSourceCountSteps(sim.integrator));
            DisplayPlot ePlot = new DisplayPlot();
//            ePlot.setUnit(new SimpleUnit(Energy.DIMENSION,46*Joule.UNIT.toSim(1)/Constants.AVOGADRO*1000,"energy","symbol",true));
            history.setDataSink(ePlot.getDataSet().makeDataSink());
            ePlot.setLabel("Energy");
            graphic.add(ePlot);
            graphic.getDisplayBox(graphic.getSimulation().getBox(0)).setPixelUnit(new Pixel(25));
            graphic.makeAndDisplayFrame();

            List<DataPump> pumps = graphic.getController().getDataStreamPumps();
            final MeterPotentialEnergy meterPE = new MeterPotentialEnergy(sim.potentialMaster);
            meterPE.setBox(sim.box);

            MeterDADBWaterTIP4P meterDADB = new MeterDADBWaterTIP4P(sim.space, meterPE, sim.potentialMaster, Kelvin.UNIT.toSim(temperature), sim.latticeCoordinates);
            meterDADB.doRotation = doRotation;
            meterDADB.doTranslation = doTranslation;
            DataFork forkPE = new DataFork();
            AccumulatorHistory historyPE = new AccumulatorHistory(new HistoryCollapsingAverage());
            DataPumpListener pumpPE = new DataPumpListener(meterPE, null, 10);
            pumps.add(pumpPE);
            sim.integrator.getEventManager().addListener(pumpPE);
            AccumulatorAverageCollapsing accPE = new AccumulatorAverageCollapsing(100, 1);
            accPE.setPushInterval(1);
            forkPE.addDataSink(accPE);
            forkPE.addDataSink(historyPE);

            DisplayTextBoxesCAE displayPE = new DisplayTextBoxesCAE();
            displayPE.setAccumulator(accPE);
            graphic.add(displayPE);

            DataFork forkHMA = new DataFork();
            AccumulatorHistory historyHMA = new AccumulatorHistory(new HistoryCollapsingAverage());
            DataPumpListener pumpHMA = new DataPumpListener(meterDADB, forkHMA, 10);
            pumps.add(pumpHMA);

            sim.integrator.getEventManager().addListener(pumpHMA);
            AccumulatorAverageCollapsing accHMA = new AccumulatorAverageCollapsing(100, 1);
            accHMA.setPushInterval(1);
            forkHMA.addDataSink(historyHMA);
            forkHMA.addDataSink(accHMA);

            DisplayTextBoxesCAE displayHMA = new DisplayTextBoxesCAE();
            displayHMA.setAccumulator(accHMA);
            graphic.add(displayHMA);

            DisplayPlot plotPE = new DisplayPlot();
            historyHMA.addDataSink(plotPE.getDataSet().makeDataSink());
            historyPE.addDataSink(plotPE.getDataSet().makeDataSink());
            plotPE.setLabel("PE");
            plotPE.setLegend(new DataTag[]{meterPE.getTag()}, "Conv");
            plotPE.setLegend(new DataTag[]{meterDADB.getTag()}, "HMA");
            graphic.add(plotPE);

            if (false) {
                MeterDADBWaterTIP4P.justU = true;
                pumpPE.setDataSink(historyPE);
            } else {
                DataProcessor processorAnh = new DataProcessor() {
                    protected DataInfoDouble dataInfo;
                    protected DataDouble data = new DataDouble();
                    protected DataTag tag = new DataTag();

                    public DataPipe getDataCaster(IDataInfo inputDataInfo) {
                        return null;
                    }

                    protected IDataInfo processDataInfo(IDataInfo inputDataInfo) {
                        dataInfo = new DataInfoDouble("conv U anh", Null.DIMENSION);
                        dataInfo.addTags(inputDataInfo.getTags());
                        dataInfo.addTag(tag);
                        return dataInfo;
                    }

                    protected IData processData(IData inputData) {
                        int N = sim.box.getMoleculeList().getMoleculeCount();
                        double fac = (doRotation ? 1.5 : 0) * N + (doTranslation ? 1.5 : 0) * (N - 1);
                        data.x = inputData.getValue(0) - latticeEnergy - fac * Kelvin.UNIT.toSim(temperature);
                        return data;
                    }
                };
                pumpPE.setDataSink(processorAnh);
                processorAnh.setDataSink(forkPE);
            }

            return;
        }//Graphic!!


        final MeterPotentialEnergyFromIntegrator meterPE = new MeterPotentialEnergyFromIntegrator(sim.integrator);

        int blockSize = numSteps >= 1000 ? (numSteps / 1000) : 1;

        MeterDADBWaterTIP4P meterDADB = new MeterDADBWaterTIP4P(sim.space, meterPE, sim.potentialMaster, Kelvin.UNIT.toSim(temperature), sim.latticeCoordinates);
//        meterDADB.getData();
//        MeterDADB.justU = true;
        meterDADB.doTranslation = doTranslation;
        meterDADB.doRotation = doRotation;

        AccumulatorAverageFixed accumulatorAverageFixedDADB = new AccumulatorAverageFixed(blockSize);
        DataPumpListener dataPumpListenerDADB = new DataPumpListener(meterDADB, accumulatorAverageFixedDADB, 4 * numCells * 46);

        AccumulatorAverageFixed accumulatorAverageFixedPE = new AccumulatorAverageFixed(blockSize);
        DataPumpListener dataPumpListenerPE = new DataPumpListener(meterPE, accumulatorAverageFixedPE, 4);


        sim.ai.setMaxSteps(numSteps / 10);
        sim.getController().actionPerformed();

        //TODO
//        meterDADB.debug();
//        System.out.println("Do debug in meterDADB");

        sim.ai.setMaxSteps(numSteps);
        sim.integrator.getEventManager().addListener(dataPumpListenerDADB);
        sim.integrator.getEventManager().addListener(dataPumpListenerPE);
        sim.getController().reset();
        sim.getController().actionPerformed();

        long endTime = System.currentTimeMillis();
        DateFormat date = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
        Calendar cal = Calendar.getInstance();
        System.out.println(date.format(cal.getTime()));
        System.out.println("Time taken (in mins): " + (endTime - startTime) / (1000.0 * 60.0));
        System.out.println("numSteps = " + numSteps);

        double mappingAverage = accumulatorAverageFixedDADB.getData(AccumulatorAverage.AVERAGE).getValue(0);
        double mappingError = accumulatorAverageFixedDADB.getData(AccumulatorAverage.ERROR).getValue(0);
        double mappingCor = accumulatorAverageFixedDADB.getData(AccumulatorAverage.BLOCK_CORRELATION).getValue(0);

        System.out.println("rCutLJ = " + rCutLJ);
        System.out.println("rCutRealES = " + rCutRealES);
        System.out.println("temperature = " + temperature);
        IMoleculeList molecules = sim.box.getMoleculeList();
        System.out.println("beginLE = " + latticeEnergy);
        int N = sim.box.getMoleculeList().getMoleculeCount();
        double fac = (doRotation ? 1.5 : 0) * N + (doTranslation ? 1.5 : 0) * (N - 1);
        System.out.println("harmonicE = " + fac * Kelvin.UNIT.toSim(temperature));
        System.out.println("main:doTranslation:   " + doTranslation + "  doRotation:  " + doRotation);

        System.out.println("mappingAverage= " + mappingAverage + "   mappingError= " + mappingError + " mappingCor= " + mappingCor);

        double PEAverage = accumulatorAverageFixedPE.getData(AccumulatorAverage.AVERAGE).getValue(0);
        double PEAError = accumulatorAverageFixedPE.getData(AccumulatorAverage.ERROR).getValue(0);
        double PECor = accumulatorAverageFixedPE.getData(AccumulatorAverage.BLOCK_CORRELATION).getValue(0);


        System.out.println("PEAverage= " + (PEAverage - latticeEnergy - fac * Kelvin.UNIT.toSim(temperature)) +
                "  PEeError= " + PEAError + "  PECor= " + PECor);

        System.out.println("PE-lattice = " + (PEAverage - latticeEnergy));
//        ConfigurationFile config = new ConfigurationFile(numCells + "ncFinalPos");
//        config.initializeCoordinates(sim.box);
//        final double endLatticeEnergy = meterPE2.getDataAsScalar();
//        System.out.println("endLE = " + endLatticeEnergy);


    }

    public static class WaterOrientationDefinition implements MoleculeSiteSource.MoleculeOrientationDefinition {
        protected final OrientationFull3D or;
        protected final Vector v1, v2;

        public WaterOrientationDefinition(Space space) {
            or = new OrientationFull3D(space);
            v1 = space.makeVector();
            v2 = space.makeVector();

        }

        public IOrientation getOrientation(IMolecule molecule) {
            IAtomList leafList = molecule.getChildList();
            Vector h1 = leafList.getAtom(0).getPosition();
            Vector h2 = leafList.getAtom(1).getPosition();
            Vector o = leafList.getAtom(2).getPosition();
            Vector m = leafList.getAtom(3).getPosition();
            v1.Ev1Mv2(m, o);
            v1.normalize();
            v2.Ev1Mv2(h2, h1);
            v2.normalize();
            or.setDirections(v1, v2);
            //v1 is a0 and v2 is a
            return or;
        }
    }

    public static class WaterDADBParam extends ParameterBase {
        public int numCells = 1;
        public int numSteps = 1000;
        public double temperature = 10;
        public double rCutLJ = 11;
        public double rCutRealES = 11;
        public double kCut = 1.5;
        public boolean isIce = false;
        public boolean runGraphic = false;
        public boolean unitCells = false;
        public boolean doRotation = true;
        public boolean doTranslation = false;
    }
}
