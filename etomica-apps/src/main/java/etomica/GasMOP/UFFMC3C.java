package etomica.GasMOP;

import etomica.action.BoxInflate;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.AccumulatorHistory;
import etomica.data.DataPumpListener;
import etomica.data.DataSourceCountSteps;
import etomica.data.history.HistoryCollapsingAverage;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.types.DataGroup;
import etomica.graphics.*;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.lattice.LatticeCubicFcc;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.potential.*;
import etomica.potential.UFF.*;
import etomica.simulation.Simulation;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesManager;
import etomica.units.Kelvin;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

import java.util.*;


public class UFFMC3C extends Simulation {
    public Box box;
    public IntegratorMC integrator;
    public SpeciesManager sm;
    ISpecies speciesMOP, autoMOP;
    public UFFMC3C() {
        this(new SimParamsUFF());
    }
    /**
     * Creates simulation with default parameters from
    /**
     * Creates simulation with the given parameters
**/
    public UFFMC3C(SimParamsUFF params) {
        super(Space3D.getInstance());
        PDBReaderMOP pdbReaderMOP = new PDBReaderMOP();
        SetPotential SetPotential = new SetPotential();


        if (params.ifMOPPresent) {
            speciesMOP = pdbReaderMOP.getSpecies(params.confName1, false, new Vector3D(0.0, 0.0, 0.0), false, false);
            addSpecies(speciesMOP);
        }else if(params.ifautoMOP){
            autoMOP = pdbReaderMOP.getSpeciesMOP(params.confName1, false, new Vector3D(0.0, 0.0, 0.0), false);
            addSpecies(autoMOP);
        }

        box = this.makeBox();
        ArrayList<ArrayList<Integer>> connectivityModified1 = pdbReaderMOP.getConnectivityModified();
        Map<Integer, String> atomMapModified1 = pdbReaderMOP.getAtomMapModified();
        double truncCal =params.rc;
        DistCalc distCalc = new DistCalc();
        if(params.ifMOPPresent){
            Map<Integer, Vector> listPositions = pdbReaderMOP.getPositions();
            box.getBoundary().setBoxSize(params.boxsize);
            box.setNMolecules(speciesMOP, params.numAtoms);
        } else {
            box.getBoundary().setBoxSize(new Vector3D(60,60,60));
        }
        //System.out.println(boxSize);
        if(params.ifautoMOP){
            distCalc.getAutoMOPBox(Space3D.getInstance(), params.struc, params.autoStruc, speciesMOP, box, connectivityModified1, atomMapModified1);
        }
        Map<Double, List<Integer[]>>distMap = new HashMap<>();
        //tetra= 6, cube = 12 octahedron = 12, dodecahedron = 30, icosahedron = 30;
        sm = new SpeciesManager.Builder().addSpecies(speciesMOP).build();

        //Bonding Parameters GasOne
        UniversalSimulation.makeAtomPotentials(sm);
        ArrayList<ArrayList<Integer>> connectedAtoms1 = pdbReaderMOP.getConnectivity();
        Map<Integer,String> atomMap1 = pdbReaderMOP.getAtomMapWithoutRunning();
        ArrayList<Integer> bondList1 = pdbReaderMOP.getBondList(connectedAtoms1, atomMap1);
        Map<String, double[]> atomicPotMap1 = pdbReaderMOP.atomicPotMap();
        ArrayList<Integer> bondsNum1 = pdbReaderMOP.getBonds();
        Map<Integer, String> atomIdentifierMapModified1 = pdbReaderMOP.getModifiedAtomIdentifierMap();
        List<int[]> dupletsSorted1= pdbReaderMOP.getDupletesSorted();
        List<int[]>tripletsSorted1= pdbReaderMOP.getAnglesSorted();
        List<int[]>quadrupletsSorted1= pdbReaderMOP.getTorsionSorted();
        Map<String[],List<int[]>> bondTypesMap1= pdbReaderMOP.idenBondTypes(dupletsSorted1, atomIdentifierMapModified1);
        Map<String[],List<int[]>> angleTypesMap1= pdbReaderMOP.idenAngleTypes(tripletsSorted1, atomIdentifierMapModified1);
        Map<String[],List<int[]>> torsionTypesMap1= pdbReaderMOP.idenTorsionTypes(quadrupletsSorted1, atomIdentifierMapModified1);
        //MOP
        PotentialMasterBonding pmBonding = new PotentialMasterBonding(sm, box);
        PotentialMasterBonding.FullBondingInfo bondingInfo1 = new PotentialMasterBonding.FullBondingInfo(sm);
       // SetPotential.setBondStretch(speciesMOP, bondTypesMap1, angleTypesMap1, torsionTypesMap1,bondsNum1,bondList1, quadrupletsSorted1, atomIdentifierMapModified1,atomicPotMap1, bondingInfo1,  pmBonding);
        SetPotential.setupBondStrech(speciesMOP, bondTypesMap1, angleTypesMap1, torsionTypesMap1,bondsNum1,bondList1, quadrupletsSorted1, atomIdentifierMapModified1,atomicPotMap1, bondingInfo1, pmBonding);
        PotentialMasterCell potentialMasterCell = new PotentialMasterCell(getSpeciesManager(), box, 5, pmBonding.getBondingInfo());

        //GasOne
        List<AtomType> listGas;
        List<List<AtomType>> listMOPGasPairs = null;

        //NonBonded
        potentialMasterCell.doAllTruncationCorrection = false;
        HashSet<AtomType> uniqueGasAtomTypes = new HashSet<AtomType>();
        uniqueGasAtomTypes.addAll(speciesMOP.getAtomTypes());
        System.out.println(uniqueGasAtomTypes);
        listMOPGasPairs = SetPotential.listFinal(uniqueGasAtomTypes);



        //MOP-Gas

        // System.out.println(listMOPGas);

        LJUFF[] p2LJMOPGas = new LJUFF[listMOPGasPairs.size()];
        P2Electrostatic[] p2ElectroMOPGas = new P2Electrostatic[listMOPGasPairs.size()];
        IPotential2[] p2mopgas = new IPotential2[listMOPGasPairs.size()];
        //   LJCOMPASS[] p2LJMOPGasCOMPASS = new LJCOMPASS[listMOPGasPairs.size()];

        if(params.ifautoMOP || params.ifMOPPresent){
            SetPotential.doLJElectrostatic(listMOPGasPairs, potentialMasterCell, p2LJMOPGas, p2ElectroMOPGas, p2mopgas, listMOPGasPairs.size(), truncCal, false);
        }

        integrator = new IntegratorMC(potentialMasterCell, random, Kelvin.UNIT.toSim(params.temperature), box);
        potentialMasterCell.doOneTruncationCorrection = true;
        potentialMasterCell.init();

        BoxInflate inflater = new BoxInflate(box, space, params.density);
        inflater.actionPerformed();
        ConfigurationLattice config = new ConfigurationLattice(new LatticeCubicFcc(space), space);
        config.initializeCoordinates(box);

        MCMoveAtom mcMoveAtom = new MCMoveAtom(random, potentialMasterCell, box);
        integrator.getMoveManager().addMCMove(mcMoveAtom);

        this.getController().addActivity(new ActivityIntegrate(integrator));
    }

    public static void main(String[] args) {

        SimParamsUFF params = new SimParamsUFF();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        } else {
            // modify parameters here for interactive testing
        }

        UFFMC3C sim = new UFFMC3C(params);
        long steps = params.steps;
        int interval = 10;
        int blocks = 100;
        long blockSize = steps / (interval * blocks);

        System.out.println("Lennard-Jones Monte Carlo simulation");
        System.out.println("N: " + params.numAtoms);
        System.out.println("T: " + params.temperature);
        System.out.println("density: " + params.density);
        System.out.println("steps: " + params.steps);

        // equilibration
        long t1 = System.currentTimeMillis();
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, steps / 10));
        System.out.println("equilibration finished");

        // data collection
        MeterPotentialEnergyFromIntegrator meterPE = new MeterPotentialEnergyFromIntegrator(sim.integrator);
        AccumulatorAverageFixed acc = new AccumulatorAverageFixed(blockSize);
        DataPumpListener pump = new DataPumpListener(meterPE, acc, interval);
        sim.integrator.getEventManager().addListener(pump);

        sim.integrator.resetStepCount();
        sim.integrator.getMoveManager().setEquilibrating(false);
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, steps));

        long t2 = System.currentTimeMillis();

        DataGroup dataPE = (DataGroup) acc.getData();
        int numAtoms = sim.getBox(0).getLeafList().size();
        double avg = dataPE.getValue(acc.AVERAGE.index) / numAtoms;
        double err = dataPE.getValue(acc.ERROR.index) / numAtoms;
        double cor = dataPE.getValue(acc.BLOCK_CORRELATION.index);

        System.out.println("energy avg: " + avg + "  err: " + err + "  cor: " + cor);
        System.out.println("time: " + (t2 - t1) * 0.001);
    }

    public static class Graphic {
        public static void main(String[] args) {
            SimParamsUFF params = new SimParamsUFF();
            if (args.length > 0) {
                ParseArgs.doParseArgs(params, args);
            } else {
                // modify parameters here for interactive testing
            }

            UFFMC3C sim = new UFFMC3C(params);
            SimulationGraphic graphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);

            DataSourceCountSteps timeSource = new DataSourceCountSteps(sim.integrator);
            MeterPotentialEnergyFromIntegrator meterPE = new MeterPotentialEnergyFromIntegrator(sim.integrator);
            AccumulatorHistory accPE = new AccumulatorHistory(new HistoryCollapsingAverage());
            accPE.setTimeDataSource(timeSource);
            DataPumpListener pumpPE = new DataPumpListener(meterPE, accPE, 10);
            sim.integrator.getEventManager().addListener(pumpPE);

            DisplayPlot historyPE = new DisplayPlot();
            accPE.setDataSink(historyPE.getDataSet().makeDataSink());
            historyPE.setLabel("PE");
            graphic.add(historyPE);

            graphic.makeAndDisplayFrame();

        }
    }

    public static class SimParamsUFF extends ParameterBase {
        public long steps = 1000000;
        public String confName1 ="F://Avagadro//molecule//co2";
        public String struc ="tetra";
        public String autoStruc = "edge";
        public boolean ifMOPPresent = true;
        public double rc = 10;
        public boolean ifautoMOP = false;
        public double density = 0.00001;
        public double temperature = 273;
        public Vector boxsize = new Vector3D(60,60,60);
        public int numAtoms = 1;
    }

}
