package etomica.models.clathrates;
import java.awt.Color;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.List;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.*;
import etomica.box.Box;
import etomica.config.ConfigurationFile;
import etomica.config.ConfigurationFileBinary;
import etomica.data.*;
import etomica.data.history.HistoryCollapsingAverage;
import etomica.data.meter.*;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.exception.ConfigurationOverlapException;
import etomica.graphics.*;
 import etomica.integrator.IntegratorVelocityVerletRattle;
import etomica.integrator.IntegratorVelocityVerletShake.BondConstraints;
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
import etomica.potential.Potential2SoftSphericalLS;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.IOrientation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.OrientationFull3D;
import etomica.space3d.Space3D;
import etomica.units.Calorie;
import etomica.units.Joule;
import etomica.units.Kelvin;
import etomica.units.Mole;
import etomica.units.Pixel;
import etomica.units.dimensions.Null;
import etomica.util.Constants;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.util.random.RandomMersenneTwister;
public class WaterDADB extends Simulation {

    public IntegratorVelocityVerletRattle integrator;
    public PotentialMaster potentialMaster;
    public SpeciesWater4P species;
    public Box box;
    public ActivityIntegrate ai;
    protected Potential2SoftSphericalLS potentialLJLS;
    protected final MoleculeAgentManager latticeCoordinates;

    public WaterDADB(final Space space, double temperature, int numCells, double rCutRealES, double rCutLJ, boolean isIce, double kCut, double shakeTol, boolean unitCells, final boolean doTranslation, final boolean doRotation, double timeInterval) {
        super(space);
//        setRandom(new RandomMersenneTwister(2));
//        System.out.println("The randomSeed is fixed");

        box = new Box(space);
        addBox(box);
        species = new SpeciesWater4P(getSpace(), true);
        addSpecies(species);
        box.setNMolecules(species, 46 * numCells * numCells * numCells);
        box.setDensity(46 / 12.03 / 12.03 / 12.03);
        ChargeAgentSourceRPM agentSource = new ChargeAgentSourceRPM(species, isIce);
        AtomLeafAgentManager<EwaldSummation.MyCharge> atomAgentManager = new AtomLeafAgentManager<EwaldSummation.MyCharge>(agentSource, box);
        int[] nC = new int[]{numCells, numCells, numCells};
        if (unitCells) {
            numCells = 1;
        }
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

        EwaldSummation potentialES = new EwaldSummation(box, atomAgentManager, space, kCut, rCutRealES);

        P2LennardJones potentialLJ = new P2LennardJones(space, sigma, epsilon);
        potentialLJLS = new Potential2SoftSphericalLS(space, rCutLJ, rC, potentialLJ);
//		potentialLJ =  new P2SoftSphericalTruncated(space, potentialLJ, rC);
        potentialMaster = new PotentialMaster();
        potentialMaster.addPotential(potentialES, new AtomType[0]);
        potentialMaster.addPotential(potentialLJLS, new AtomType[]{species.getOxygenType(), species.getOxygenType()});

        //ADD MINIMIZATION TO IT LATER

        int maxIterations = 100;
        integrator = new IntegratorVelocityVerletRattle(this, potentialMaster, box);
        integrator.setShakeTolerance(shakeTol);
        double lOH = ConformationWaterTIP4P.bondLengthOH;
        double lHH = Math.sqrt(2 * lOH * lOH * (1 - Math.cos(ConformationWaterTIP4P.angleHOH)));
        double lOM = ConformationWaterTIP4P.rOM;
//        double lMH = Math.sqrt(lOH * lOH + lOM * lOM - 2 * lOH * lOM * Math.cos(0.5 * ConformationWaterTIP4P.angleHOH));
        BondConstraints bondConstraints = new BondConstraints(new int[][]{{0, 2}, {1, 2}, {0, 1}}, new double[]{lOH, lOH, lHH}) {

            Vector vectorSum = space.makeVector();
            Vector oVector = space.makeVector();
            Vector centerMass = space.makeVector();
            Vector h1Vector = space.makeVector();
            Vector h2Vector = space.makeVector();
            Vector mVector = space.makeVector();
            //            Vector newForce = space.makeVector();
//            Vector newTorque = space.makeVector();
            Vector dr = space.makeVector();

            public void redistributeForces(IMolecule molecule, AtomLeafAgentManager<Vector> agentManager) {
                IAtomList leafList = molecule.getChildList();
                Vector h1 = leafList.get(0).getPosition();
                Vector h2 = leafList.get(1).getPosition();
                Vector o = leafList.get(2).getPosition();
                Vector m = leafList.get(3).getPosition();
                Vector h1torque = space.makeVector();
                Vector h2torque = space.makeVector();
                Vector h1h2 = space.makeVector();
                Vector om = space.makeVector();
                Vector oTorque = space.makeVector();
                Vector mTorque = space.makeVector();
                Vector totalTorque = space.makeVector();
                Vector totalForce = space.makeVector();
                double hMass = leafList.get(0).getType().getMass();
                double oMass = leafList.get(2).getType().getMass();
                double hMassPercent = hMass / (2 * hMass + oMass);
                double oMassPercent = oMass / (2 * hMass + oMass);
                centerMass.Ea1Tv1(hMass, h1);
                centerMass.PEa1Tv1(hMass, h2);
                centerMass.PEa1Tv1(oMass, o);
                centerMass.TE(1 / (2 * hMass + oMass));
                mVector.Ev1Mv2(m, centerMass);
                h1Vector.Ev1Mv2(h1, centerMass);
                h2Vector.Ev1Mv2(h2, centerMass);
                oVector.Ev1Mv2(o, centerMass);
                vectorSum.E(h1Vector);
                vectorSum.PE(h2Vector);
                vectorSum.PEa1Tv1(-2.0, oVector);
                Vector h1Force = agentManager.getAgent(leafList.get(0));
                Vector h2Force = agentManager.getAgent(leafList.get(1));
                Vector oForce = agentManager.getAgent(leafList.get(2));
                Vector mForce = agentManager.getAgent(leafList.get(3));

                boolean testForceAndTorque = false;
                Vector totalTorqueNew = space.makeVector();
                Vector totalForceNew = space.makeVector();
                if (testForceAndTorque) {
                    totalForce.E(h1Force);
                    totalForce.PE(h2Force);
                    totalForce.PE(mForce);
                    totalForce.PE(oForce);// test for total force

                    h1torque.E(h1Force);
                    h1torque.XE(h1Vector);
                    h2torque.E(h2Force);
                    h2torque.XE(h2Vector);
                    oTorque.E(oForce);
                    oTorque.XE(oVector);
                    mTorque.E(mForce);
                    mTorque.XE(mVector);
                    totalTorque.E(h1torque);
                    totalTorque.PE(h2torque);
                    totalTorque.PE(oTorque);
                    totalTorque.PE(mTorque);//test for total torque
                }

//#####################WHAT IS HAPPENING BELOW?????
                h1h2.Ev1Mv2(h2, h1);
                om.Ev1Mv2(m, o);
                Vector a0 = space.makeVector();
                Vector a1 = space.makeVector();
                Vector a2 = space.makeVector();
                a0.E(mVector);
                a0.normalize();
                a1.Ea1Tv1(-1, h1h2);   //#########WHAT DOES Ea1Tv1 DO?????
                a1.normalize();
                a2.E(a0);//equal
                a2.XE(a1);//cross a2=a2Xa1

                //separate mForce into 3 direction
                Vector f0 = space.makeVector();
                Vector f1 = space.makeVector();
                Vector f2 = space.makeVector();
                f0.Ea1Tv1(mForce.dot(a0), a0);
                f1.Ea1Tv1(mForce.dot(a1), a1);
                f2.Ea1Tv1(mForce.dot(a2), a2);

                //Translation part
                h1Force.PEa1Tv1(hMassPercent, mForce);
                h2Force.PEa1Tv1(hMassPercent, mForce);
                oForce.PEa1Tv1(oMassPercent, mForce);

                //Rotation part
                dr.Ea1Tv1(f1.dot(a1) * mVector.dot(a0) / (-2 * oVector.dot(a0) + 2 * h1Vector.dot(a0)), a1);
//                System.out.println("part1="+ f1.dot(mVector)/(-2*oVector.dot(a0)+2*h1Vector.dot(a0)));
//                System.out.println("part1Force={"+dr.getX(0)+","+dr.getX(1)+","+dr.getX(2)+"};" );
                h1Force.PE(dr);
                h2Force.PE(dr);
                oForce.PEa1Tv1(-2, dr);

                dr.Ea1Tv1(f2.dot(a2) * mVector.dot(a0) / (-2 * oVector.dot(a0) + 2 * h1Vector.dot(a0)), a2);
//                System.out.println("part2="+f2.dot(mVector)/(-2*oVector.dot(a0)+2*h1Vector.dot(a0)));
//                System.out.println("part2Force={"+dr.getX(0)+","+dr.getX(1)+","+dr.getX(2)+"};" );
                h1Force.PE(dr);
                h2Force.PE(dr);
                oForce.PEa1Tv1(-2, dr);

                mForce.E(0);//remove mForce at the end

                if (testForceAndTorque) {
                    //test for total force
                    totalForceNew.E(h1Force);
                    totalForceNew.PE(h2Force);
                    totalForceNew.PE(mForce);
                    totalForceNew.PE(oForce);//New total force

                    System.out.println("totalForceOld ={ " + totalForce.getX(0) + "," + totalForce.getX(1) + "," + totalForce.getX(2) + "};");
                    System.out.println("totalForceNew ={ " + totalForceNew.getX(0) + "," + totalForceNew.getX(1) + "," + totalForceNew.getX(2) + "};");


                    //test for total torque
                    h1torque.E(h1Force);
                    h1torque.XE(h1Vector);
                    h2torque.E(h2Force);
                    h2torque.XE(h2Vector);
                    oTorque.E(oForce);
                    oTorque.XE(oVector);
                    mTorque.E(mForce);
                    mTorque.XE(mVector);
                    totalTorqueNew.E(h1torque);
                    totalTorqueNew.PE(h2torque);
                    totalTorqueNew.PE(mTorque);
                    totalTorqueNew.PE(oTorque);//New total torque
                    System.out.println("totalTorqueOld ={ " + totalTorque.getX(0) + "," + totalTorque.getX(1) + "," + totalTorque.getX(2) + "};");
                    System.out.println("totalTorqueNew ={ " + totalTorqueNew.getX(0) + "," + totalTorqueNew.getX(1) + "," + totalTorqueNew.getX(2) + "};");
                    System.exit(2);
                }


                //redistribute the forces s.t. only translation
                if (!doRotation) {
                    totalForce.E(h1Force);
                    totalForce.PE(h2Force);
                    totalForce.PE(mForce);
                    totalForce.PE(oForce);
                    h1Force.E(0);
                    h2Force.E(0);
                    oForce.E(0);
                    mForce.E(0);
                    h1Force.PEa1Tv1(hMassPercent, totalForce);
                    h2Force.PEa1Tv1(hMassPercent, totalForce);
                    oForce.PEa1Tv1(oMassPercent, totalForce);
                }


                //redistribute s.t. only rotation
                if (!doTranslation) {
                    totalForce.E(h1Force);
                    totalForce.PE(h2Force);
                    totalForce.PE(mForce);
                    totalForce.PE(oForce);
                    h1Force.PEa1Tv1(-1.0 * hMassPercent, totalForce);
                    h2Force.PEa1Tv1(-1.0 * hMassPercent, totalForce);
                    oForce.PEa1Tv1(-1.0 * oMassPercent, totalForce);
                }

            }

            public void relaxMolecule(IMolecule molecule) {
                IAtomList leafList = molecule.getChildList();

                Vector h1 = leafList.get(0).getPosition();
                Vector h2 = leafList.get(1).getPosition();
                Vector o = leafList.get(2).getPosition();
                Vector m = leafList.get(3).getPosition();
                m.E(h1);
                m.PE(h2);
                m.PEa1Tv1(-2, o);
                m.TE(0.15 / Math.sqrt(m.squared())); // TIP4P
//            	m.TE(0.1577/Math.sqrt(m.squared())); TIP4P/Ice
                m.PE(o);

            }

        };
        integrator.getShakeAgentManager().setAgent(species, bondConstraints);
        integrator.setTimeStep(timeInterval);
        integrator.setMaxIterations(maxIterations);
//        integrator.setOrientAtom((IAtom)((IMolecule)box.getMoleculeList(speciesOrient).getAtom(0)).getChildList().getAtom(0));
        integrator.setIsothermal(true);
        integrator.setTemperature(Kelvin.UNIT.toSim(temperature));
        integrator.setThermostatInterval(100);
        integrator.setThermostatNoDrift(true);

        try {
            integrator.reset();
        } catch (ConfigurationOverlapException e) {
        }
        ai = new ActivityIntegrate(integrator);
//        System.out.println("using rigid with dt="+dt);
        getController().addAction(ai);
//        System.out.println("h1 at "+((IAtomPositioned)box.getLeafList().getAtom(0)).getPosition());
//        System.out.println("o at "+((IAtomPositioned)box.getLeafList().getAtom(2)).getPosition());

    }

    public static void main(String[] args) {
//    	System.out.println(Kelvin.UNIT.toSim(1));
//    	System.out.println(Joule.UNIT.toSim(1)/Constants.AVOGADRO);
//    	System.exit(2);
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
        double shakeTol = waterDADBParam.shakeTol;
        boolean uniteCells = waterDADBParam.unitCells;
        boolean doTranslation = waterDADBParam.doTranslation;
        boolean doRotation = waterDADBParam.doRotation;
        boolean runGraphic = waterDADBParam.runGraphic;
        double timeInterval = waterDADBParam.timeInterval;
        boolean doMapping = waterDADBParam.doMapping;
        boolean doConventional = waterDADBParam.doConventional;
        final WaterDADB sim = new WaterDADB(Space3D.getInstance(), temperature, numCells, rCutRealES, rCutLJ, isIce, kCut, shakeTol, uniteCells, doTranslation, doRotation, timeInterval);
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
//			sim.ai.setSleepPeriod(2);
            SimulationGraphic graphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, "Rattle", 1);
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
            MeterEnergy meterE = new MeterEnergy(sim.potentialMaster, sim.box);
            AccumulatorHistory history = new AccumulatorHistory(new HistoryCollapsingAverage());
            history.setTimeDataSource(new DataSourceCountTime(sim.integrator));
            DataProcessor processor = new DataProcessor() {
                protected DataInfoDouble dataInfo;
                protected DataDouble data = new DataDouble();
                protected DataTag tag = new DataTag();

                public DataPipe getDataCaster(IDataInfo inputDataInfo) {
                    return null;
                }

                protected IDataInfo processDataInfo(IDataInfo inputDataInfo) {
                    dataInfo = new DataInfoDouble("foo", Null.DIMENSION);
                    dataInfo.addTags(inputDataInfo.getTags());
                    dataInfo.addTag(tag);
                    return dataInfo;
                }

                protected IData processData(IData inputData) {
                    data.x = inputData.getValue(0) / temperature;
                    return data;
                }
            };
            DataPumpListener pump = new DataPumpListener(meterE, processor,10);
            processor.setDataSink(history);
            DisplayPlot ePlot = new DisplayPlot();
//            ePlot.setUnit(new SimpleUnit(Energy.DIMENSION,46*Joule.UNIT.toSim(1)/Constants.AVOGADRO*1000,"energy","symbol",true));
            history.setDataSink(ePlot.getDataSet().makeDataSink());
            sim.integrator.getEventManager().addListener(pump);
            ePlot.setLabel("Energy");
            graphic.add(ePlot);
            graphic.getDisplayBox(graphic.getSimulation().getBox(0)).setPixelUnit(new Pixel(25));
            graphic.makeAndDisplayFrame();
            System.out.println(meterE.getDataAsScalar() / (46 * Joule.UNIT.toSim(1) / Constants.AVOGADRO * 1000));
            System.out.println(meterE.getDataAsScalar());

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
  //          plotPE.setLegend(new DataTag[]{meterDADB.getTag()}, "HMA");
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
                        int N = sim.box.getMoleculeList().size();
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

        int blockSize = numSteps >= 1000 ? (numSteps / 1000) : 1;

        final MeterPotentialEnergy meterPE = new MeterPotentialEnergy(sim.potentialMaster);
        meterPE.setBox(sim.box);

        MeterDADBWaterTIP4P meterDADB ;
        AccumulatorAverageFixed accumulatorAverageFixedDADB = null;
        DataPumpListener dataPumpListenerDADB = null;
        if (doMapping) {
            meterDADB = new MeterDADBWaterTIP4P(sim.space, meterPE, sim.potentialMaster, Kelvin.UNIT.toSim(temperature), sim.latticeCoordinates);
            meterDADB.doTranslation = doTranslation;
            meterDADB.doRotation = doRotation;
            accumulatorAverageFixedDADB = new AccumulatorAverageFixed(blockSize);
            dataPumpListenerDADB = new DataPumpListener(meterDADB, accumulatorAverageFixedDADB, 10);
            //        MeterDADB.justU = true;
        }
        //TODO try lower interval 5-10

        AccumulatorAverageFixed accumulatorAverageFixedPE = null;
        DataPumpListener dataPumpListenerPE = null;
        if (doConventional) {
            accumulatorAverageFixedPE = new AccumulatorAverageFixed(blockSize);
            dataPumpListenerPE = new DataPumpListener(meterPE, accumulatorAverageFixedPE, 10);
        }

        sim.ai.setMaxSteps(numSteps / 10);
        sim.getController().actionPerformed();

//        meterDADB.debug();
//        System.out.println("Do debug in meterDADB");
//        System.exit(2);


        sim.ai.setMaxSteps(numSteps);
        if (doMapping) sim.integrator.getEventManager().addListener(dataPumpListenerDADB);
        if (doConventional) sim.integrator.getEventManager().addListener(dataPumpListenerPE);
        sim.getController().reset();
        sim.getController().actionPerformed();

        System.out.println("numSteps= " + numSteps);
        System.out.println("temperature= " + temperature);
//        System.out.println("rCutLJ = " + rCutLJ);
//        System.out.println("rCutRealES = " + rCutRealES);
//        System.out.println("shakeTol= " + shakeTol);
        IMoleculeList molecules = sim.box.getMoleculeList();
        System.out.println("beginLE = " + latticeEnergy);
        int N = sim.box.getMoleculeList().size();
        double fac = (doRotation ? 1.5 : 0) * N + (doTranslation ? 1.5 : 0) * (N - 1);
        System.out.println("harmonicE= " + fac * Kelvin.UNIT.toSim((molecules.size() * temperature)));
        System.out.println("main:doTranslation:   " + doTranslation + "  doRotation:  " + doRotation);
        System.out.println("timeInterval= " + timeInterval);


        long endTime = System.currentTimeMillis();
        DateFormat date = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
        Calendar cal = Calendar.getInstance();
        System.out.println(date.format(cal.getTime()));
        double totalTime = (endTime - startTime) / (1000.0 * 60.0);
        System.out.println("totalTime= " + totalTime + " mins");

        if (doMapping) {
            double mappingAverage = accumulatorAverageFixedDADB.getData(AccumulatorAverage.AVERAGE).getValue(0);
            double mappingError = accumulatorAverageFixedDADB.getData(AccumulatorAverage.ERROR).getValue(0);
            double mappingCor = accumulatorAverageFixedDADB.getData(AccumulatorAverage.BLOCK_CORRELATION).getValue(0);
            System.out.println("mappingAverage=\t" + mappingAverage +
                    "   mappingError=\t" + mappingError + " mappingCor=\t" + mappingCor + " time= " + totalTime);
        }

        if (doConventional) {
            double PEAverage = accumulatorAverageFixedPE.getData(AccumulatorAverage.AVERAGE).getValue(0);
            double PEAError = accumulatorAverageFixedPE.getData(AccumulatorAverage.ERROR).getValue(0);
            double PECor = accumulatorAverageFixedPE.getData(AccumulatorAverage.BLOCK_CORRELATION).getValue(0);
            System.out.println("PEAverage=\t" + (PEAverage - latticeEnergy - fac * Kelvin.UNIT.toSim(temperature)) +
                    "  PEeError=\t" + PEAError + "  PECor=\t" + PECor + " time= " + totalTime);
        }
//        System.out.println("PE-lattice= " + (PEAverage - latticeEnergy));
//        ConfigurationFile config = new ConfigurationFile(numCells + "ncFinalPos");
//        config.initializeCoordinates(sim.box);
//        final double endLatticeEnergy = meterPE2.getDataAsScalar();
//        System.out.println("endLE = " + endLatticeEnergy);



    }

    public static class WaterOrientationDefinition implements etomica.normalmode.MoleculeSiteSource.MoleculeOrientationDefinition {
        protected final OrientationFull3D or;
        protected final Vector v1, v2;

        public WaterOrientationDefinition(Space space) {
            or = new OrientationFull3D(space);
            v1 = space.makeVector();
            v2 = space.makeVector();

        }

        public IOrientation getOrientation(IMolecule molecule) {
            IAtomList leafList = molecule.getChildList();
            Vector h1 = leafList.get(0).getPosition();
            Vector h2 = leafList.get(1).getPosition();
            Vector o = leafList.get(2).getPosition();
            Vector m = leafList.get(3).getPosition();
            v1.Ev1Mv2(m, o);
            v1.normalize();
            v2.Ev1Mv2(h2, h1);
            v2.normalize();
            or.setDirections(v1, v2);
            //v1 is a0 and v2 is a1
            return or;
        }
    }

    public static class WaterDADBParam extends ParameterBase {
        public int numCells = 1;
        public int numSteps = 1000;
        public double timeInterval = 0.002;
        public double temperature = 50;
        public double rCutLJ = 11;
        public double rCutRealES = 11;
        public double kCut = 1.5;
        public boolean isIce = false;
        public double shakeTol = 1e-12;
        public boolean runGraphic = false;
        public boolean unitCells = false;
        public boolean doRotation = true;
        public boolean doTranslation = true;
        public boolean doMapping = true;
        public boolean doConventional = true;
    }

}
