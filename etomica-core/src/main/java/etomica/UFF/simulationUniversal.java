package etomica.UFF;

import etomica.action.BoxImposePbc;
import etomica.action.BoxInflate;
import etomica.action.IAction;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.DiameterHashByType;
import etomica.atom.SpeciesMethane;
import etomica.box.Box;
import etomica.config.ConfigurationFile;
import etomica.config.ConfigurationLattice;
import etomica.data.*;
import etomica.data.history.HistoryCollapsingAverage;
import etomica.data.history.HistoryCollapsingDiscard;
import etomica.data.meter.MeterDensity;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.meter.MeterPressure;
import etomica.data.types.DataDouble;
import etomica.graphics.*;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.integrator.mcmove.MCMoveMolecule;
import etomica.integrator.mcmove.MCMoveMoleculeRotate;
import etomica.integrator.mcmove.MCMoveVolume;
import etomica.lattice.LatticeCubicFcc;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.potential.*;
import etomica.potential.UFF.P3BondAngleUFF;
import etomica.potential.UFF.P4BondTorsionUFF;
import etomica.potential.UFF.UFF;
import etomica.potential.compute.PotentialCompute;
import etomica.potential.compute.PotentialComputeAggregate;
import etomica.simulation.Simulation;
import etomica.simulation.prototypes.MCMoveWiggle;
import etomica.simulation.prototypes.MeterTorsionAngle;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.SpBuilder;
import etomica.species.SpeciesGeneral;
import etomica.species.SpeciesManager;
import etomica.units.*;
import etomica.units.dimensions.Null;
import etomica.util.Constants;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.util.random.RandomMersenneTwister;

import java.awt.*;
import java.util.List;
import java.util.*;


public class simulationUniversal extends Simulation {

    public PotentialCompute pcAgg;
    public IntegratorMC integratorMC;
    public SpeciesGeneral species;
    public Box box;
    public MCMoveVolume mcMoveVolume;
    public MCMoveMolecule mcMoveMolecule;
    // static ArrayList<ArrayList<Integer>> connectivity = new ArrayList<>();
    //Map<Integer, String> atomMap = new HashMap<>();
    public simulationUniversal(Space space, double density, int numMoleules, double temperature, String configFileName, double rc, double pressure, double rcElectric) {
        super(space);
        int atom1 = 0, atom2 = 0, atom3=0;
        String confName = "F:/methane";
        String atomName1 = null, atomName2= null, atomName3= null;
        setRandom(new RandomMersenneTwister(1));
        species = (SpeciesGeneral) SpeciesMethane.buildMethane(true, confName );
        addSpecies(species);
        box = this.makeBox();
        box.setNMolecules(species, numMoleules);
        new BoxInflate(box, space, density).actionPerformed();
        SpeciesManager sm = new SpeciesManager.Builder().addSpecies(species).build();
        PotentialMasterBonding pmBonding = new PotentialMasterBonding(sm, box);

        ArrayList<ArrayList<Integer>> connectedAtoms = PDBBuilder.getConnectivityWithSpecies(confName);
        System.out.println(connectedAtoms+ ": connectedAtom");
        ArrayList<ArrayList<Integer>> connectivityModified = PDBBuilder.getconnectivityModified(connectedAtoms);
        System.out.println(connectivityModified+ ": connectedAtomModified" );
        Map<Integer,String> atomMap = PDBBuilder.getAtomMap(connectedAtoms);
        System.out.println(atomMap + ": atomMap");
        HashMap<Integer, String> atomMapModified = PDBBuilder.getatomMapModified(atomMap);
        System.out.println(atomMapModified + ": atomMapModified");
        List<int[]> duplets = PDBBuilder.getBondList(connectivityModified);
        System.out.println(Arrays.deepToString(duplets.toArray())+ ": listOfBonds");
        List<int[]> triplets = PDBBuilder.getAngleList(connectivityModified);
        System.out.println(Arrays.deepToString(triplets.toArray())+ ": listOfAngleModified");
        List<int[]> quadruplets = PDBBuilder.getTorsionList(connectedAtoms);
        System.out.println(Arrays.deepToString(quadruplets.toArray())+ " listOfTorsionModified");
        Set<String> uniqueElements = PDBBuilder.uniqueElementIdentifier();
        System.out.println(uniqueElements + "Set of Unique Elements");
        ArrayList<Integer> bondList = PDBBuilder.getBonding(confName);
        Unit kcals = new UnitRatio(new PrefixedUnit(Prefix.KILO,Calorie.UNIT),Mole.UNIT);
        Unit Radian = new Radian();
        Map<String, double[]> atomicPotMap = PDBBuilder.atomicPotMap();
        System.out.println(atomicPotMap + "atomicPotMap");
        Map<Integer, String> atomIdentifierMapModified = PDBBuilder.atomIdentifierMapModified(connectivityModified, atomMapModified);
        System.out.println(atomIdentifierMapModified + "atomIdentifierMapModified");
        Set<String> uniqueAtoms = PDBBuilder.uniqueElementIdentifier();
        System.out.println(uniqueAtoms);

        AtomType typeC_3 = species.getTypeByName("C");
        //AtomType typeC_2 = species.getTypeByName("C");
        AtomType typeH = species.getTypeByName("H");
        //AtomType typeCl = species.getTypeByName("CL");
        SpBuilder spBuilder = new SpBuilder();

        String keyC = "C-ene";
        String keyH = "H";
        String keyCl = "Cl";
        String keyO ="O";

        double[] Carb = SpBuilder.atomicPot(keyC);
        //System.out.println(Arrays.toString(Carb) + "Carb");
        double[] Hydr = SpBuilder.atomicPot(keyH);
        //System.out.println(Arrays.toString(Hydr) + "Hydr");
        double[] Chlor = SpBuilder.atomicPot(keyCl);
        //System.out.println(Arrays.toString(Chlor) + "Chlor");
        UFF uff = new UFF();
        //System.out.println("After duplet");
        List<int[]> arraytoList = new ArrayList<>();

        P2Harmonic bondUFFArray = null;
        P3BondAngleUFF angleUFF = null;

        List<String[]> bondListWithnames = new ArrayList<>();
        bondListWithnames.add(new String[]{atomName1, atomName2});

        List<int[]> bondListAtoms = new ArrayList<>();
        bondListAtoms.add(new int[]{atom1, atom2});
        // System.out.println(bondListAtoms + " " + bondListWithnames);
        //Combine both loops
        System.out.println("\n\nbonds");
        System.out.println("Entered bonding"); int i =0;
        for (int[] bond : duplets) {
            atom1 = bond[0];
            atom2 = bond[1];
            atomName1 = atomIdentifierMapModified.get(atom1);
            atomName2 = atomIdentifierMapModified.get(atom2);
            double[] atomOnePot = atomicPotMap.get(atomName1);
            double[] atomTwoPot = atomicPotMap.get(atomName2);
            bondUFFArray =  UFF.bondUFF(atomOnePot[0], atomTwoPot[0],1, 0.1332, atomOnePot[5], atomTwoPot[5], atomOnePot[6], atomTwoPot[6]);
            arraytoList.add(bond);
            i++;
        }
        //System.out.println("out of bonding");
        System.out.println(Arrays.deepToString(arraytoList.toArray()));
        pmBonding.setBondingPotentialPair(species, bondUFFArray, arraytoList);

        //System.out.println("\n \n Angles start");
        List<String[]> angleListWithNames = new ArrayList<>();
       // System.out.println("In angles");
        for (int[] angle : triplets) {
            Unit degree = Degree.UNIT;
            atom1 = angle[0];
            atom2 = angle[1];
            atom3 = angle[2];
            atomName1 = atomIdentifierMapModified.get(atom1);
            atomName2 = atomIdentifierMapModified.get(atom2);
            atomName3 = atomIdentifierMapModified.get(atom3);
            double[] atomOnePot = atomicPotMap.get(atomName1);
            double[] atomTwoPot = atomicPotMap.get(atomName2);
            double[] atomThreePot = atomicPotMap.get(atomName3);
            double thetha0Rad = degree.toSim(atomTwoPot[1]);
            arraytoList.add(angle);
            angleUFF =  UFF.angleUFF(atomOnePot[0], atomTwoPot[0], atomThreePot[0],1, 0.1332, atomOnePot[5], atomTwoPot[5], atomThreePot[5], atomOnePot[6], atomTwoPot[6],atomThreePot[6], thetha0Rad);

        }
        //System.out.println("out of angles");
        pmBonding.setBondingPotentialTriplet(species, angleUFF, triplets);

        P4BondTorsionUFF torsionUFF = null;
        double Vi =0, Vj =0, V=0, Vtrue=0,  type;
        int p;
        for(int[] quad: quadruplets){
            type = 0;
            p = 0;
            //System.out.println(Arrays.toString(quad) + "array is printed");

                    //A-B-C-D
                    String atomName = atomMapModified.get(quad[1]);
                    Vi = UFF.switchCaseTorsion(atomName);

                    //System.out.println(atomName + " "+ 1 + " "+ Vi);
                    int bondListValue = bondList.get(quad[1]);
                    //System.out.println(bondListValue);
                    p = p + bondListValue;

                   atomName = atomMapModified.get(quad[2]);
                    Vj = UFF.switchCaseTorsion(atomName);
                   // System.out.println(atomName + " "+ 2+ " "+ Vj);
                    bondListValue = bondList.get(quad[2]);
                    p = p + bondListValue;

                    V = Math.sqrt(Vi*Vj);
                    Vtrue = kcals.toSim(V);
                  //  System.out.println(V + " V " + p + " p value");
                   
                    torsionUFF = UFF.torsionUFF(Vtrue, p);
        }
        pmBonding.setBondingPotentialTriplet(species, angleUFF, quadruplets);
        
        PotentialMasterCell potentialMaster = new PotentialMasterCell(getSpeciesManager(), box, 2, pmBonding.getBondingInfo());
        potentialMaster.doAllTruncationCorrection = true;
        pcAgg = new PotentialComputeAggregate(pmBonding, potentialMaster);
        integratorMC = new IntegratorMC(pcAgg, random, temperature, box);
        getController().addActivity(new ActivityIntegrate(integratorMC));

        mcMoveMolecule = new MCMoveMolecule(random, pcAgg, box);
        integratorMC.getMoveManager().addMCMove(mcMoveMolecule);

        MCMoveMoleculeRotate rotateMove = new MCMoveMoleculeRotate(random, pcAgg, box);
        integratorMC.getMoveManager().addMCMove(rotateMove);

        MCMoveWiggle wiggleMove = new MCMoveWiggle(random, pcAgg, box);
        integratorMC.getMoveManager().addMCMove(wiggleMove);

        MCMoveAtom moveAtom = new MCMoveAtom(random, pcAgg, box);
        integratorMC.getMoveManager().addMCMove(moveAtom);


        //P4BondTorsionUFF torsionUFF = null;
        P2Electrostatic electrostaticFF = null;

        double epsilonH = kcals.toSim(Hydr[2]);
        double epsilonC = kcals.toSim(Carb[2]);;
        double sigmaH = Hydr[3];
        double sigmaC = Carb[3];
        double sigmaCH = (Hydr[3]+Carb[3])/2;
        TruncationFactory tf = new TruncationFactoryForceShift(rc);
        P2LennardJones p2LJHH = uff.vdw(epsilonH, epsilonH, sigmaH, sigmaH);
        P2LennardJones p2LJCH = uff.vdw(epsilonH, epsilonC, sigmaH, sigmaC);
        P2LennardJones p2LJCC = uff.vdw(epsilonC, epsilonC, sigmaC, sigmaC);
        IPotential2 p2HH = tf.make(p2LJHH);
        IPotential2 p2HC = tf.make(p2LJCH);
        IPotential2 p2CC = tf.make(p2LJCC);

        double chargeH = Hydr[5];
        double chargeC = Carb[5];
        TruncationFactory tfe = new TruncationFactoryForceShift(rcElectric);
        P2Electrostatic p2ElectrostaticHH = uff.electroUFF(chargeH, chargeH);
        P2Electrostatic p2ElectrostaticCH = uff.electroUFF(chargeC, chargeH);
        P2Electrostatic p2ElectrostaticCC = uff.electroUFF(chargeC, chargeC);

        IPotential2 p2EHH = tfe.make(p2ElectrostaticHH);
        p2EHH = new P2SoftSphericalSumTruncated(rc, p2EHH, p2LJHH);
        IPotential2 p2ECH = tfe.make(p2ElectrostaticCH);
        p2ECH = new P2SoftSphericalSumTruncated(rc, p2ECH, p2LJCH);
        IPotential2 p2ECC = tfe.make(p2ElectrostaticCC);
        p2ECC = new P2SoftSphericalSumTruncated(rc, p2ECC, p2LJCC);
        potentialMaster.setPairPotential(typeH, typeH, p2EHH);
        potentialMaster.setPairPotential(typeC_3, typeH, p2ECH);
        potentialMaster.setPairPotential(typeC_3, typeC_3, p2ECC);

        if (configFileName != null) {
            ConfigurationFile config = new ConfigurationFile(configFileName);
            config.initializeCoordinates(box);
            BoxImposePbc.imposePBC(box);
        }
        else {
            ConfigurationLattice configuration = new ConfigurationLattice(new LatticeCubicFcc(space), space);
            configuration.initializeCoordinates(box);
            potentialMaster.init();
            double u0 = potentialMaster.computeAll(false);
            int x = 1;
            while (u0 > 1e6*numMoleules) {
                x *= 0.99;
                p2LJCC.setSigma(x*sigmaC);
                p2LJCH.setSigma(x*sigmaCH);
                p2LJHH.setSigma(x*sigmaH);
                ((P2SoftSphericalSumTruncatedForceShifted)p2HH).setTruncationRadius(rc);
                ((P2SoftSphericalSumTruncatedForceShifted)p2HC).setTruncationRadius(rc);
                ((P2SoftSphericalSumTruncatedForceShifted)p2CC).setTruncationRadius(rc);
                u0 = potentialMaster.computeAll(false);
                //System.out.println("Inside looop");
            }
            integratorMC.reset();
            while (u0 > 1e4*numMoleules) {
                while (u0 > 1e4 * numMoleules) {
                    integratorMC.doStep();
                    u0 = integratorMC.getPotentialEnergy();
                }
                while (x < 1 && u0 <= 1e4 * numMoleules) {
                    x /= 0.99;
                    if (x > 1) x = 1;
                    p2LJCC.setSigma(x * sigmaC);
                    p2LJCH.setSigma(x * sigmaCH);
                    p2LJHH.setSigma(x * sigmaH);
                    ((P2SoftSphericalSumTruncatedForceShifted) p2HH).setTruncationRadius(rc);
                    ((P2SoftSphericalSumTruncatedForceShifted) p2HC).setTruncationRadius(rc);
                    ((P2SoftSphericalSumTruncatedForceShifted) p2CC).setTruncationRadius(rc);
                    u0 = potentialMaster.computeAll(false);
                }
                integratorMC.reset();
            }
        }
    }

    public static void main(String[] args) {
        final String APP_NAME = "Methane Universal";
        OctaneParams params = new OctaneParams();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        }
        else {
            params.numSteps = 3500000;
            params.density = 0.00003;
            params.configFilename = null; // "octane";
            params.graphics = true;
        }


        Unit dUnit = new SimpleUnit(Null.DIMENSION, 1/(72.15/ Constants.AVOGADRO*1e24), "Density", "g/cm^3", false);

        double temperatureK = params.temperatureK;
        double temperature = Kelvin.UNIT.toSim(temperatureK);
        System.out.println("Tsim "+temperature);
        int numMolecules = params.numMolecules;
        double density = dUnit.toSim(params.density);
        boolean graphics = params.graphics;
        long numSteps = params.numSteps;
        String configFilename = params.configFilename;
        double rc = params.rc;
        double rcElectric = params.rcElectric;
        double pressureKPa = params.pressureKPa;
        Unit pUnit = new PrefixedUnit(Prefix.KILO, Pascal.UNIT);
        double pressure = pUnit.toSim(pressureKPa);

        System.out.println(numSteps+" steps");
        System.out.println("rc: "+rc);
        System.out.println("pressure "+ pressureKPa);
        System.out.println("initial density "+ density);
        System.out.println("initial density (g/cm^3) "+ dUnit.fromSim(density));

        final simulationUniversal sim = new simulationUniversal(Space3D.getInstance(), density, numMolecules, temperature, configFilename, rc, pressure, rcElectric);

        MeterPotentialEnergyFromIntegrator meterU = new MeterPotentialEnergyFromIntegrator(sim.integratorMC);
        sim.integratorMC.getPotentialCompute().init();
        sim.integratorMC.reset();
        System.out.println("u0/N "+(meterU.getDataAsScalar()/numMolecules));

        MeterPressure meterP = new MeterPressure(sim.box, sim.pcAgg);
        meterP.setTemperature(temperature);
        meterP.doCallComputeAll(true);
        DataProcessorForked dpZ = new DataProcessorForked() {
            final DataDouble.DataInfoDouble dataInfo = new DataDouble.DataInfoDouble("Z", Null.DIMENSION);
            final DataDouble data = new DataDouble();

            @Override
            protected IData processData(IData inputData) {
                data.x = inputData.getValue(0) / temperature / density;
                return data;
            }

            @Override
            protected IDataInfo processDataInfo(IDataInfo inputDataInfo) {
                return dataInfo;
            }
        };
        //System.out.println("Reached before dataforked");
        DataProcessorForked dpZm1oR = new DataProcessorForked() {
            DataDouble.DataInfoDouble dataInfo = new DataDouble.DataInfoDouble("(Z-1)/rho", Null.DIMENSION);
            DataDouble data = new DataDouble();

            @Override
            protected IData processData(IData inputData) {
                data.x = (inputData.getValue(0) / temperature / density - 1) / density;
                return data;
            }

            @Override
            protected IDataInfo processDataInfo(IDataInfo inputDataInfo) {
                return dataInfo;
            }
        };
        DataFork forkP = new DataFork(new IDataSink[]{dpZ, dpZm1oR});

        if (true) {
            sim.getController().addActivity(new ActivityIntegrate(sim.integratorMC, 5000000));
            sim.getController().addActionSequential(new IAction() {
                @Override
                public void actionPerformed() {
                    sim.integratorMC.getMoveManager().addMCMove(sim.mcMoveVolume);
                }
            });
            sim.getController().addActivity(new ActivityIntegrate(sim.integratorMC));
            final SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, "Universal MC", 3);
            //System.out.println("Reached after simulation graphic");
            DiameterHashByType dhbt = (DiameterHashByType) simGraphic.getDisplayBox(sim.box).getDiameterHash();
          /*  dhbt.setDiameter(sim.species.getAtomType(0), 3.75);
            dhbt.setDiameter(sim.species.getAtomType(1), 3.95);*/
            ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.species.getTypeByName("C"), Color.WHITE);

            ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.species.getTypeByName("H"), Color.BLUE);
            //System.out.println("Simulation starts before Datasource");
            DataSourceCountSteps timeSource = new DataSourceCountSteps(sim.integratorMC);

            simGraphic.getController().getReinitButton().setPostAction(simGraphic.getPaintAction(sim.box));

            List<DataPump> dataPumps = simGraphic.getController().getDataStreamPumps();
           // System.out.println("Reached dataPump");
            DataSourceCountSteps timer = new DataSourceCountSteps(sim.integratorMC);
            DisplayTextBox timerBox = new DisplayTextBox();
            timerBox.setLabel("Steps");
            DataPumpListener pumpSteps = new DataPumpListener(timer, timerBox, numMolecules);
            sim.integratorMC.getEventManager().addListener(pumpSteps);
            simGraphic.add(timerBox);

            MeterDensity meterDensity = new MeterDensity(sim.box());
            AccumulatorHistory accDensity = new AccumulatorHistory(new HistoryCollapsingAverage());
            accDensity.setTimeDataSource(timer);
            DataPumpListener pumpDensity = new DataPumpListener(meterDensity, accDensity, 10);
            sim.integratorMC.getEventManager().addListener(pumpDensity);
            dataPumps.add(pumpDensity);

            DisplayPlot historyDensity = new DisplayPlot();
            accDensity.setDataSink(historyDensity.getDataSet().makeDataSink());
            historyDensity.setLabel("Density");
            historyDensity.setUnit(dUnit);
            simGraphic.add(historyDensity);

            Unit perN = new SimpleUnit(Null.DIMENSION, numMolecules, "1/N", "1/N", false);

            AccumulatorHistory historyU = new AccumulatorHistory(new HistoryCollapsingDiscard());
            historyU.setTimeDataSource(timer);
            AccumulatorHistory historyU2 = new AccumulatorHistory(new HistoryCollapsingAverage());
            historyU2.setTimeDataSource(timer);
            AccumulatorAverageCollapsing avgEnergy = new AccumulatorAverageCollapsing();
            avgEnergy.setPushInterval(10);
            DataFork forkU = new DataFork(new IDataSink[]{historyU, historyU2, avgEnergy});
            DataPumpListener pumpU = new DataPumpListener(meterU, forkU, 10);
            dataPumps.add(pumpU);
            sim.integratorMC.getEventManager().addListener(pumpU);
            DisplayPlotXChart plotU = new DisplayPlotXChart();
            plotU.setLabel("U");
            historyU.addDataSink(plotU.makeSink("U"));
            plotU.setLegend(new DataTag[]{historyU.getTag()}, "samples");
            historyU2.addDataSink(plotU.makeSink("Uavg"));
            plotU.setLegend(new DataTag[]{historyU2.getTag()}, "avg");
            plotU.setUnit(perN);
            simGraphic.add(plotU);

            simGraphic.getController().getDataStreamPumps().add(pumpU);

            DisplayTextBoxesCAE display = new DisplayTextBoxesCAE();
            display.setAccumulator(avgEnergy);
            display.setUnit(perN);
            simGraphic.add(display);

            meterP.setTemperature(temperature);
            meterP.doCallComputeAll(true);
            AccumulatorHistory historyP = new AccumulatorHistory(new HistoryCollapsingDiscard());
            historyP.setTimeDataSource(timer);
            AccumulatorHistory historyP2 = new AccumulatorHistory(new HistoryCollapsingAverage());
            historyP2.setTimeDataSource(timer);
            AccumulatorAverageCollapsing avgP = new AccumulatorAverageCollapsing();
            forkP.addDataSink(historyP);
            forkP.addDataSink(historyP2);
            forkP.addDataSink(avgP);
            DataPumpListener pumpP = new DataPumpListener(meterP, forkP, numMolecules);
            dataPumps.add(pumpP);
            sim.integratorMC.getEventManager().addListener(pumpP);
            DisplayPlotXChart plotP = new DisplayPlotXChart();
            plotP.setLabel("P");
            historyP.addDataSink(plotP.makeSink("P"));
            plotP.setLegend(new DataTag[]{historyP.getTag()}, "samples");
            historyP2.addDataSink(plotP.makeSink("Pavg"));
            plotP.setLegend(new DataTag[]{historyP2.getTag()}, "avg");
            plotP.setUnit(pUnit);
            simGraphic.add(plotP);
            simGraphic.getController().getDataStreamPumps().add(pumpP);

            DisplayTextBoxesCAE displayP = new DisplayTextBoxesCAE();
            displayP.setAccumulator(avgP);
            simGraphic.add(displayP);

            AccumulatorHistory historyZ = new AccumulatorHistory(new HistoryCollapsingDiscard());
            historyZ.setTimeDataSource(timer);
            AccumulatorHistory historyZ2 = new AccumulatorHistory(new HistoryCollapsingAverage());
            historyZ2.setTimeDataSource(timer);
            AccumulatorAverageCollapsing avgZ = new AccumulatorAverageCollapsing();
            dpZ.addDataSink(historyZ);
            dpZ.addDataSink(historyZ2);
            dpZ.addDataSink(avgZ);
            DisplayPlotXChart plotZ = new DisplayPlotXChart();
            plotZ.setLabel("Z");
            historyZ.addDataSink(plotZ.makeSink("Zsamples"));
            plotZ.setLegend(new DataTag[]{historyZ.getTag()}, "samples");
            historyZ2.addDataSink(plotZ.makeSink("Zavg"));
            plotZ.setLegend(new DataTag[]{historyZ2.getTag()}, "avg");
            simGraphic.add(plotZ);

            DisplayTextBoxesCAE displayZ = new DisplayTextBoxesCAE();
            displayZ.setLabel("Z");
            displayZ.setAccumulator(avgZ);
            simGraphic.add(displayZ);

            AccumulatorHistory historyZ_ = new AccumulatorHistory(new HistoryCollapsingDiscard());
            historyZ_.setTimeDataSource(timer);
            AccumulatorHistory historyZ_2 = new AccumulatorHistory(new HistoryCollapsingAverage());
            historyZ_2.setTimeDataSource(timer);
            AccumulatorAverageCollapsing avgZ_ = new AccumulatorAverageCollapsing();
            dpZm1oR.addDataSink(historyZ_);
            dpZm1oR.addDataSink(historyZ_2);
            dpZm1oR.addDataSink(avgZ_);
            DisplayPlotXChart plotZ_ = new DisplayPlotXChart();
            plotZ_.setLabel("(Z-1)/rho");
            historyZ_.addDataSink(plotZ_.makeSink("samples"));
            plotZ_.setLegend(new DataTag[]{historyZ_.getTag()}, "samples");
            historyZ_2.addDataSink(plotZ_.makeSink("avg"));
            plotZ_.setLegend(new DataTag[]{historyZ_2.getTag()}, "avg");
            simGraphic.add(plotZ_);

            DisplayTextBoxesCAE displayZ_ = new DisplayTextBoxesCAE();
            displayZ_.setLabel("(Z-1)/rho");
            displayZ_.setAccumulator(avgZ_);
            simGraphic.add(displayZ_);
            simGraphic.makeAndDisplayFrame(APP_NAME);


            MeterTorsionAngle meterTorsion = new MeterTorsionAngle(sim.box, 2, 3, 4, 5);
            AccumulatorAverageBlockless accTorsion1 = new AccumulatorAverageBlockless();
            AccumulatorAverageCollapsing accTorsion2 = new AccumulatorAverageCollapsing();
            AccumulatorHistory historyTorsion = new AccumulatorHistory(new HistoryCollapsingDiscard());
            historyTorsion.setTimeDataSource(timer);
            AccumulatorHistory historyTorsion2 = new AccumulatorHistory(new HistoryCollapsingAverage());
            historyTorsion2.setTimeDataSource(timer);
            DataFork forkTorsion = new DataFork(new IDataSink[]{accTorsion1, accTorsion2, historyTorsion, historyTorsion2});
            DataPumpListener pumpTorsion = new DataPumpListener(meterTorsion, forkTorsion, numMolecules);
            dataPumps.add(pumpTorsion);
            sim.integratorMC.getEventManager().addListener(pumpTorsion);
            DisplayTextBoxesCAE displayTorsion = new DisplayTextBoxesCAE();
            displayTorsion.setLabel("Torsion cos");
            displayTorsion.setAccumulator(accTorsion2);
            DisplayPlotXChart plotTorsion = new DisplayPlotXChart();
            plotTorsion.setLabel("torsion");
            historyTorsion.addDataSink(plotTorsion.makeSink("samples"));
            plotTorsion.setLegend(new DataTag[]{historyTorsion.getTag()}, "samples");
            historyTorsion2.addDataSink(plotTorsion.makeSink("avg"));
            plotTorsion.setLegend(new DataTag[]{historyTorsion2.getTag()}, "avg");
            simGraphic.add(plotTorsion);
            simGraphic.makeAndDisplayFrame();
            return;



        }

        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integratorMC, numSteps/100));
        // sim.integratorMC.getMoveManager().addMCMove(sim.mcMoveMolecule);

        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integratorMC, numSteps/5));

        long samples = numSteps / (numMolecules* 8L);
        long bs = samples / 100;
        if (bs == 0) bs = 1;

        AccumulatorAverageFixed accU = new AccumulatorAverageFixed((numSteps/10)/100);
        DataPumpListener pumpU = new DataPumpListener(meterU, accU, 10);
        sim.integratorMC.getEventManager().addListener(pumpU);

        AccumulatorAverageFixed accP = new AccumulatorAverageFixed(bs);
        forkP.addDataSink(accP);
        AccumulatorAverageFixed accZ = new AccumulatorAverageFixed(bs);
        dpZ.addDataSink(accZ);
        AccumulatorAverageFixed accZm1oR = new AccumulatorAverageFixed(bs);
        dpZm1oR.addDataSink(accZm1oR);
        DataPumpListener pumpP = new DataPumpListener(meterP, forkP, 8*numMolecules);
        sim.integratorMC.getEventManager().addListener(pumpP);

        long t1 = System.nanoTime();
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integratorMC, numSteps));
        long t2 = System.nanoTime();

        IData dataU = accU.getData();
        double avgU = dataU.getValue(AccumulatorAverage.AVERAGE.index) / numMolecules;
        double errU = dataU.getValue(AccumulatorAverage.ERROR.index) / numMolecules;
        double corU = dataU.getValue(AccumulatorAverage.BLOCK_CORRELATION.index);
        System.out.println("U: " +" "+avgU+ "   err: "+" "+errU+"   cor: " +" "+corU);

        IData dataP = accP.getData();
        double avgP = dataP.getValue(AccumulatorAverage.AVERAGE.index);
        double errP = dataP.getValue(AccumulatorAverage.ERROR.index);
        double corP = dataP.getValue(AccumulatorAverage.BLOCK_CORRELATION.index);
        System.out.println("P: " +" "+avgP+ "   err: " +" "+errP+ "   cor: " +corP);

        IData dataZ = accZ.getData();
        double avgZ = dataZ.getValue(AccumulatorAverage.AVERAGE.index);
        double errZ = dataZ.getValue(AccumulatorAverage.ERROR.index);
        double corZ = dataZ.getValue(AccumulatorAverage.BLOCK_CORRELATION.index);
        System.out.println("Z: "+" "+avgZ+"   err: "+" "+errZ+"   cor: "+" "+corZ);

        IData dataZ_ = accZm1oR.getData();
        double avgZ_ = dataZ_.getValue(AccumulatorAverage.AVERAGE.index);
        double errZ_ = dataZ_.getValue(AccumulatorAverage.ERROR.index);
        double corZ_ = dataZ_.getValue(AccumulatorAverage.BLOCK_CORRELATION.index);
        System.out.println("(Z-1)/rho: "+" "+avgZ_+"   err: "+" "+errZ_+"   cor: "+" "+corZ_);

        System.out.println("time: "+" "+(t2-t1)/1e9);
        /*
        0.06147
         */

    }

    public static class OctaneParams extends ParameterBase {
        public double temperatureK = 300;
        public int numMolecules = 100;
        public double density = 0.00005;
        public boolean graphics = false;
        public long numSteps = 200000;
        public String configFilename = null;
        public double rc = 14;
        public double rcElectric = 20;
        public double pressureKPa = 1402;
    }
}

