package etomica.GasMOP;

import etomica.action.BoxInflate;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.DiameterHashByType;
import etomica.box.Box;
import etomica.box.RandomPositionSource;
import etomica.box.RandomPositionSourceRectangular;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.DataPumpListener;
import etomica.data.histogram.HistogramVectorEnergy;
import etomica.data.histogram.HistogramVectorSimple;
import etomica.data.meter.*;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.Integrator;
import etomica.integrator.IntegratorEvent;
import etomica.integrator.IntegratorListener;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.*;
import etomica.math.DoubleRange;
import etomica.molecule.CenterOfMass;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.molecule.MoleculeSourceRandomMolecule;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.potential.*;
import etomica.potential.TraPPE.SpeciesGasTraPPE;
import etomica.potential.UFF.*;
import etomica.potential.compute.PotentialComputeAggregate;
import etomica.potential.compute.PotentialComputeEwaldFourier;
import etomica.simulation.Simulation;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesManager;
import etomica.units.*;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.util.collections.IntArrayList;
import etomica.util.random.IRandom;
import etomica.util.random.RandomMersenneTwister;
import etomica.virial.mcmove.MCMoveClusterAngle;
import etomica.virial.mcmove.MCMoveClusterAngleBend;

import java.awt.*;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.List;

import static etomica.potential.TraPPE.SpeciesGasTraPPE.ChemForm.*;

public class GCMCGO extends Simulation {
    public IntegratorMC integrator;
    public MCMoveMolecule mcMoveMolecule;
    public PotentialMasterCell potentialMasterCell;
    public MCMoveMoleculeRotate mcMoveMoleculeRotate;
    //  public MCMoveInsertDelete mcMoveID, mcMoveIDTwo;
    public ISpecies speciesMOP, speciesGas, speciesGrapheneOne, speciesGrapheneFive,speciesSecondGas, autoMOP ;
    public ISpecies speciesMOPTwo, speciesMOPThree,speciesMOPFour, speciesMOPFive, speciesCOF;
    public Box box, boxNew;
    public P2LennardJones potential;
    public SpeciesManager sm;
    public static int atom1, atom2, atom3;
    public int numCarbon;
    public static double sigmaIKey, sigmaJKey, epsilonIKey, epsilonJKey;
    public int i;
    public static String atomName1,atomName2, atomName3 ;
    public static List<List<AtomType>> listMOPGasMixed = new ArrayList<>();
    public static List<List<AtomType>> listGrapheneGasMixed = new ArrayList<>();
    public static double massMembrane;
    public static HistogramVectorSimple targHist;
    public GCMCGO(String confNameGasOne,String confNameGraphene,Vector grapheneOne, Vector grapheneTwo, Vector grapheneThree, Vector grapheneFour, int numMolOne, double temperature, double truncatedRadius, double truncatedRadiusLJ, double sigma, double mu1, double mu2, boolean ifGraphenePresent, boolean ifMultipleGraphenePresent, boolean ifMoveRotateMoves, Vector boxSize, boolean makeAllMove, boolean doElectrostatics, boolean ifGrapheneDatafile, boolean isGasCOMPASS, boolean isGasTraPPE, boolean isGasOPLS) throws IOException {
        super(Space3D.getInstance());
        //  System.out.println(struc +" "+ autoStruc);
        //Make Species

        PDBReaderMOP pdbReaderMOP = new PDBReaderMOP();
        PDBReaderMOP pdbReaderMOP2 = new PDBReaderMOP();
        PDBReaderMOP pdbReaderCOF = new PDBReaderMOP();
        AutoMOP autoMOP1 = new AutoMOP();
        PDBReaderReplica pdbReaderReplica = new PDBReaderReplica();
        PDBReaderReplica pdbReaderReplicaNew = new PDBReaderReplica();
        GeneralGrapheneReader grapheneReader = new GeneralGrapheneReader();
        GrapheneReaderXYZPDB grapheneReaderXYZPDB = new GrapheneReaderXYZPDB();
        GasOPLS gasOPLS = new GasOPLS();
        SetPotential SetPotential = new SetPotential();
        GCMCGOParams gcmcgoParams = new GCMCGOParams();

        ArrayList<ArrayList<Integer>> connectivity1 = new ArrayList<>();
        Map<Integer,String> atomMap1 = new HashMap<>();

        double  truncCal = truncatedRadius;
        GetSpeciesMOP getSpeciesMOP = new GetSpeciesMOP();


        if (ifGraphenePresent){
            if(ifGrapheneDatafile){
                speciesGrapheneOne = grapheneReaderXYZPDB.getSpecies(confNameGraphene,new Vector3D(0, 0, 0), false);
                addSpecies(speciesGrapheneOne);
            }else {
                speciesGrapheneOne = grapheneReader.getSpecies(confNameGraphene, new Vector3D(0, 0, 0), false);
                addSpecies(speciesGrapheneOne);
            }
            if (isGasTraPPE){
                SpeciesGasTraPPE speciesGasTraPPE = new SpeciesGasTraPPE();
                if (confNameGasOne.equals("ch4")) {
                    SpeciesGasTraPPE.ChemForm = CH4;
                } else if (confNameGasOne.equals("ethane")) {
                    SpeciesGasTraPPE.ChemForm = C2H6;
                } else if (confNameGasOne.equals("propane")) {
                    SpeciesGasTraPPE.ChemForm = C3H8;
                } else if (confNameGasOne.equals("ethene")) {
                    SpeciesGasTraPPE.ChemForm = C2H4;
                } else if (confNameGasOne.equals("propene")) {
                    SpeciesGasTraPPE.ChemForm = C3H6;
                } else if (confNameGasOne.equals("co2")) {
                    SpeciesGasTraPPE.ChemForm = CO2;
                } else if (confNameGasOne.equals("o2")) {
                    SpeciesGasTraPPE.ChemForm = O2;
                } else if (confNameGasOne.equals("n2")) {
                    SpeciesGasTraPPE.ChemForm = N2;
                }else if (confNameGasOne.equals("nh3")) {
                    SpeciesGasTraPPE.ChemForm = NH3;
                }else if (confNameGasOne.equals("1butene")) {
                    SpeciesGasTraPPE.ChemForm = C4H8ene_1;
                }else if (confNameGasOne.equals("13butadiene")) {
                    SpeciesGasTraPPE.ChemForm = C4H8ene_13;
                }else if (confNameGasOne.equals("nbutane")) {
                    SpeciesGasTraPPE.ChemForm = C4H10;
                }else if (confNameGasOne.equals("methylpropene")) {
                    SpeciesGasTraPPE.ChemForm = methylPropene_2;
                }else if (confNameGasOne.equals("methylpropane")) {
                    SpeciesGasTraPPE.ChemForm = methylPropane_2;
                }else if (confNameGasOne.equals("transbutene")) {
                    SpeciesGasTraPPE.ChemForm = transButene;
                }else if (confNameGasOne.equals("cisbutene")) {
                    SpeciesGasTraPPE.ChemForm = cisButene;
                }
                speciesGas = speciesGasTraPPE.speciesGasTraPPE(Space3D.getInstance(), SpeciesGasTraPPE.ChemForm, false);
                addSpecies(speciesGas);
            } else if (isGasOPLS) {
                speciesGas = gasOPLS.getSpecies(confNameGasOne, new Vector3D(0,0,0),false );
                addSpecies(speciesGas);
            }
        }else {
            if (isGasTraPPE){
                SpeciesGasTraPPE speciesGasTraPPE = new SpeciesGasTraPPE();
                if (confNameGasOne.equals("ch4")) {
                    SpeciesGasTraPPE.ChemForm = CH4;
                } else if (confNameGasOne.equals("ethane")) {
                    SpeciesGasTraPPE.ChemForm = C2H6;
                } else if (confNameGasOne.equals("propane")) {
                    SpeciesGasTraPPE.ChemForm = C3H8;
                } else if (confNameGasOne.equals("ethene")) {
                    SpeciesGasTraPPE.ChemForm = C2H4;
                } else if (confNameGasOne.equals("propene")) {
                    SpeciesGasTraPPE.ChemForm = C3H6;
                } else if (confNameGasOne.equals("co2")) {
                    SpeciesGasTraPPE.ChemForm = CO2;
                } else if (confNameGasOne.equals("o2")) {
                    SpeciesGasTraPPE.ChemForm = O2;
                } else if (confNameGasOne.equals("n2")) {
                    SpeciesGasTraPPE.ChemForm = N2;
                }else if (confNameGasOne.equals("nh3")) {
                    SpeciesGasTraPPE.ChemForm = NH3;
                }else if (confNameGasOne.equals("1butene")) {
                    SpeciesGasTraPPE.ChemForm = C4H8ene_1;
                }else if (confNameGasOne.equals("13butadiene")) {
                    SpeciesGasTraPPE.ChemForm = C4H8ene_13;
                }else if (confNameGasOne.equals("nbutane")) {
                    SpeciesGasTraPPE.ChemForm = C4H10;
                }else if (confNameGasOne.equals("methylpropene")) {
                    SpeciesGasTraPPE.ChemForm = methylPropene_2;
                }else if (confNameGasOne.equals("methylpropane")) {
                    SpeciesGasTraPPE.ChemForm = methylPropane_2;
                }else if (confNameGasOne.equals("transbutene")) {
                    SpeciesGasTraPPE.ChemForm = transButene;
                }else if (confNameGasOne.equals("cisbutene")) {
                    SpeciesGasTraPPE.ChemForm = cisButene;
                }
                speciesGas = speciesGasTraPPE.speciesGasTraPPE(Space3D.getInstance(), SpeciesGasTraPPE.ChemForm, false);
                addSpecies(speciesGas);
            } else if (isGasOPLS) {
                speciesGas = gasOPLS.getSpecies(confNameGasOne, new Vector3D(0,0,0),false );
                addSpecies(speciesGas);
            }
        }
        box = this.makeBox();
        box.getBoundary().setBoxSize(boxSize);

        Map<Double, List<Integer[]> >distMap = new HashMap<>();
        //tetra= 6, cube = 12 octahedron = 12, dodecahedron = 30, icosahedron = 30;
        if(ifGraphenePresent ){
            if(ifMultipleGraphenePresent){
                List<Vector> oldPositions = new ArrayList<>();
                massMembrane = speciesGrapheneOne.getMass()*numMolOne;
                box.setNMolecules(speciesGrapheneOne, 2);
                IMolecule moleculeMOPZero = box.getMoleculeList().get(0);
                while (oldPositions.size() < moleculeMOPZero.getChildList().size()) {
                    oldPositions.add(space.makeVector());
                }
                Vector vecgrapheneOne = new Vector3D(0.0,0.0, -10);
                Vector finalVecgrapheneOne = vecgrapheneOne;
                System.out.println(finalVecgrapheneOne);
                moleculeMOPZero.getChildList().forEach(atom -> {
                    oldPositions.get(atom.getIndex()).E(atom.getPosition());
                    atom.getPosition().PE(finalVecgrapheneOne);
                    Vector shift = box.getBoundary().centralImage(atom.getPosition());
                    atom.getPosition().PE(shift);
                });
                vecgrapheneOne = new Vector3D(0.0,0.0, 10);
                System.out.println(vecgrapheneOne);
                IMolecule moleculeMOPOne = box.getMoleculeList().get(1);
                Vector finalVecgrapheneOne1 = vecgrapheneOne;
                moleculeMOPOne.getChildList().forEach(atom -> {
                    oldPositions.get(atom.getIndex()).E(atom.getPosition());
                    atom.getPosition().PE(finalVecgrapheneOne1);
                    Vector shift = box.getBoundary().centralImage(atom.getPosition());
                    atom.getPosition().PE(shift);
                });
            } else {
                box.addNewMolecule(speciesGrapheneOne);
            }
        }
        Map<Integer, String> atomIdentifierMapModified1 = pdbReaderMOP.getModifiedAtomIdentifierMap();
        if (ifGraphenePresent){
            sm = new SpeciesManager.Builder().addSpecies(speciesGrapheneOne).addSpecies(speciesGas).build();
        } else {
            sm = new SpeciesManager.Builder().addSpecies(speciesGas).build();
        }

        ArrayList<ArrayList<Integer>> connectedAtoms2= null;
        ArrayList<ArrayList<Integer>> connectivityModified2= null;
        Map<Integer,String> atomMap2= null;
        ArrayList<Integer> bondList2= null;
        Map<String, double[]> atomicPotMap2= null;
        Map<Integer, String> atomIdentifierMapModified2= null;
        List<int[]>dupletsSorted2= null;
        List<int[]>tripletsSorted2= null;
        List<int[]>quadrupletsSorted2= null;
        ArrayList<Integer> bondsNum2 = null;
        Map<String[], List<int[]>> bondTypesMap2 = null;
        Map<String[], List<int[]>> angleTypesMap2 = null;
        Map<String[], List<int[]>> torsionTypesMap2 = null;
        int[] harmonicBond = null;
        List<int[]> bonds = new ArrayList<>();
        List<int[]> angles = new ArrayList<>();
        List<int[]> torsions = new ArrayList<>();
        //Bonding Parameters GasOne
        if (isGasCOMPASS && !isGasTraPPE) {
            connectedAtoms2 = pdbReaderReplica.getConnectivityWithoutRunning();
            connectivityModified2 = pdbReaderReplica.getConnectivityModifiedWithoutRunning();
            atomMap2 = pdbReaderReplica.getAtomMapWithoutRunning();
            bondList2 = pdbReaderReplica.getBondList(connectedAtoms2, atomMap2);
            atomicPotMap2 = pdbReaderReplica.atomicPotMap();
            atomIdentifierMapModified2 = pdbReaderReplica.getatomIdentifierMapModified();
            dupletsSorted2 = pdbReaderReplica.getDupletesSorted();
            tripletsSorted2 = pdbReaderReplica.getAnglesSorted();
            quadrupletsSorted2 = pdbReaderReplica.getTorsionSorted();
            bondsNum2 = pdbReaderReplica.getBonds();
            bondTypesMap2= pdbReaderReplica.idenBondTypes(dupletsSorted2, atomIdentifierMapModified2);
            angleTypesMap2= pdbReaderReplica.idenAngleTypes(tripletsSorted2, atomIdentifierMapModified2);
            torsionTypesMap2= pdbReaderReplica.idenTorsionTypes(quadrupletsSorted2, atomIdentifierMapModified2);
        } else {
            connectedAtoms2 = pdbReaderReplica.getConnectivity();
            connectivityModified2 = pdbReaderReplica.getConnectivityModifiedWithoutRunning();
            atomMap2 = pdbReaderReplica.getAtomMapWithoutRunning();
            bondList2 = pdbReaderReplica.getBondList(connectedAtoms2, atomMap2);
            atomicPotMap2 = pdbReaderReplica.atomicPotMap();
            atomIdentifierMapModified2 = pdbReaderReplica.getatomIdentifierMapModified();
            dupletsSorted2 = pdbReaderReplica.getDupletesSorted();
            tripletsSorted2 = pdbReaderReplica.getAnglesSorted();
            quadrupletsSorted2 = pdbReaderReplica.getTorsionSorted();
            bondsNum2 = pdbReaderReplica.getBonds();
            bondTypesMap2 = pdbReaderReplica.idenBondTypes(dupletsSorted2, atomIdentifierMapModified2);
            angleTypesMap2 = pdbReaderReplica.idenAngleTypes(tripletsSorted2, atomIdentifierMapModified2);
            torsionTypesMap2 = pdbReaderReplica.idenTorsionTypes(quadrupletsSorted2, atomIdentifierMapModified2);
        }

        //SetBonding Potential
        UniversalSimulation.makeAtomPotentials(sm);
        PotentialMasterBonding pmBonding = new PotentialMasterBonding(sm, box);
        PotentialMasterCell potentialMasterCell = new PotentialMasterCell(getSpeciesManager(), box, 5, pmBonding.getBondingInfo());

       // PotentialComputeEwaldFourier ewaldFourier = new PotentialComputeEwaldFourier(getSpeciesManager(), box);
       // PotentialComputeEwaldFourier.EwaldParams params = ewaldFourier.getOptimalParams(2.62, 0);
      //  System.out.println(params);
        //GasOne
        Map<Integer, double[]> anglePotential = new HashMap<>();
        Map<Integer, double[]> bondPotential = new HashMap<>();
        Map<Integer, double[]> torsionPotential = new HashMap<>();
        GasOPLS gasOPLS1 = new GasOPLS();
        if (!isGasTraPPE && !isGasOPLS) {
            SetPotential.setBondStretch(speciesGas, bondTypesMap2, angleTypesMap2, torsionTypesMap2, bondsNum2, bondList2, quadrupletsSorted2, atomIdentifierMapModified2, atomicPotMap2, pmBonding);
        } else if (isGasOPLS) {
            String[] parts = confNameGasOne.split("//");
            String gasName = parts[parts.length-1];
            gasOPLS1.makePotential(gasName,  pmBonding);
        } else {
            String[] parts = confNameGasOne.split("//");
            String gasName = parts[parts.length-1];
            //SetPotential.setBondStretchTraPPE(speciesGas, bonds, pmBonding, gasName);
            SetPotential.setBondStretchTraPPE(speciesGas, bonds, angles, torsions, pmBonding, gasName);
        }

        List<AtomType> listGas;
        //NonBonded
        potentialMasterCell.doAllTruncationCorrection = false;

        listGas = speciesGas.getUniqueAtomTypes();
        Map<String, Double> chargeCoeffGas = new HashMap<>();
        Map<String, double[]> potentialCoeffGas = new HashMap<>();
        Map<String, Double> chargeCoeffGraph = new HashMap<>();
        Map<String, double[]> potentialCoeffGraph = new HashMap<>();
        if (isGasOPLS && ifGrapheneDatafile){
            chargeCoeffGas = gasOPLS.getChargeCoeff();
            chargeCoeffGraph = grapheneReaderXYZPDB.getChargeCoeff();
            potentialCoeffGas = gasOPLS.getCoeffPotential();
            potentialCoeffGraph = grapheneReaderXYZPDB.getCoeffPotential();
        }
        //Graphene-Gas
        if (ifGraphenePresent) {
            List<List<AtomType>> listGrapheneGasFinal = SetPotential.listGrapheneSpecial(speciesGrapheneOne.getUniqueAtomTypes(), listGas);
            //  List<List<AtomType>> listGrapheneGasFinal = SetPotential.listFinal(listGrapheneGasUnique);
            P2Electrostatic[] p2ElectroGasgraphene = new P2Electrostatic[listGrapheneGasFinal.size()];
            IPotential2[] p2electroGasgraphene = new IPotential2[listGrapheneGasFinal.size()];
            LJUFF[] p2LJGrapheneGas = new LJUFF[listGrapheneGasFinal.size()];
            if (isGasTraPPE) {
                SetPotential.doLJElectrostatic(listGrapheneGasFinal, potentialMasterCell, p2LJGrapheneGas, p2ElectroGasgraphene, listGrapheneGasFinal.size(), truncatedRadius, false, true);
            } else if (isGasOPLS) {
                SetPotential.setPotentialGrapheneGas(listGrapheneGasFinal, potentialMasterCell, p2LJGrapheneGas, p2ElectroGasgraphene, listGrapheneGasFinal.size(), truncatedRadius, doElectrostatics, isGasTraPPE, isGasOPLS, ifGrapheneDatafile, chargeCoeffGas, potentialCoeffGas, chargeCoeffGraph, potentialCoeffGraph );
            }else {
                SetPotential.doLJElectrostatic(listGrapheneGasFinal, potentialMasterCell, p2LJGrapheneGas, p2ElectroGasgraphene, p2electroGasgraphene, listGrapheneGasFinal.size(), truncatedRadiusLJ, doElectrostatics);
            }
        }


        List<List<AtomType>> listGasGasPairs = null;
        listGas = speciesGas.getUniqueAtomTypes();
        listGasGasPairs = SetPotential.listFinal(listGas);
        //NonBonded

        P2Electrostatic[] p2ElectroGasGas = new P2Electrostatic[listGasGasPairs.size()];
        IPotential2[] p2electroGasGas = new IPotential2[listGasGasPairs.size()];
        LJUFF[] p2LJGasGas = new LJUFF[listGasGasPairs.size()];
        if (isGasTraPPE || isGasOPLS){
            SetPotential.setPotentialGasGasGO(listGasGasPairs, potentialMasterCell, p2LJGasGas, p2ElectroGasGas, listGasGasPairs.size(), truncatedRadius, doElectrostatics, isGasTraPPE, isGasOPLS, ifGrapheneDatafile, potentialCoeffGas, chargeCoeffGas);
        }else {
            SetPotential.doLJElectrostatic(listGasGasPairs, potentialMasterCell, p2LJGasGas, p2ElectroGasGas, p2electroGasGas, listGasGasPairs.size(), truncatedRadiusLJ, doElectrostatics);
        }


        potentialMasterCell.doOneTruncationCorrection = true;
        potentialMasterCell.init();
        PotentialComputeAggregate pcAgg = new PotentialComputeAggregate(potentialMasterCell, pmBonding);
        integrator = new IntegratorMC(pcAgg, random, Kelvin.UNIT.toSim(temperature), box);
        if (ifMoveRotateMoves) {
            mcMoveMolecule = new MCMoveMolecule(random, potentialMasterCell, box);
            mcMoveMolecule.setStepSize( sigma);
            ((MCMoveStepTracker) mcMoveMolecule.getTracker()).setTunable(false);
            mcMoveMoleculeRotate = new MCMoveMoleculeRotate(random, potentialMasterCell, box);
            mcMoveMoleculeRotate.setStepSize( sigma);
            ((MCMoveStepTracker) mcMoveMoleculeRotate.getTracker()).setTunable(false);
            //((MoleculeSourceRandomMolecule) mcMoveMolecule.getMoleculeSource()).setSpecies(speciesGas);
            ((MoleculeSourceRandomMolecule) mcMoveMoleculeRotate.getMoleculeSource()).setSpecies(speciesGas);
            integrator.getMoveManager().addMCMove(mcMoveMoleculeRotate);
            if (makeAllMove) {
                ((MoleculeSourceRandomMolecule) mcMoveMolecule.getMoleculeSource()).setSpecies(speciesMOP);
                ((MoleculeSourceRandomMolecule) mcMoveMoleculeRotate.getMoleculeSource()).setSpecies(speciesMOP);
            }
            ((MoleculeSourceRandomMolecule) mcMoveMolecule.getMoleculeSource()).setSpecies(speciesGas);
            ((MoleculeSourceRandomMolecule) mcMoveMoleculeRotate.getMoleculeSource()).setSpecies(speciesGas);
        }

        MCMoveInsertDelete mcMoveID = new MCMoveInsertDelete(potentialMasterCell, random, space);
        mcMoveID.setBox(box);
        mcMoveID.setMu(mu1);
        mcMoveID.setSpecies(speciesGas);
        integrator.getMoveManager().addMCMove(mcMoveID);
    }

    public static void main(String[] args) throws IOException {
        GCMCGOParams params = new GCMCGOParams();
        ParseArgs.doParseArgs(params, args);
        int numAtomOne = params.numAtomOne;
        int numAtomTwo = params.numAtomTwo;
        boolean ifautoMOP = params.ifautoMOP;
        String confNameGraphene = params.confNameGraphene;
        Vector boxSize = params.boxSize;
        double truncatedRadius = params.truncatedRadius;
        double truncatedRadiusLJ =0.0;
        double temperature = params.temperature;
        Vector centreMOP = params.centreMOP;
        Vector grapheneOne = params.grapheneOne;
        Vector grapheneTwo = params.grapheneTwo;
        Vector grapheneThree = params.grapheneThree;
        Vector grapheneFour = params.grapheneFour;
        Vector centreMOPTwo = params.centreMOPTwo;
        boolean ifGraphenePresent = params.ifGraphenePresent;
        boolean ifSecondGasPresent = params.ifSecondGasPresent;
        boolean ifMultipleGraphenePresent = params.ifMultipleGraphenePresent;
        boolean ifMoveRotateMoves = params.ifMoveRotateMoves;
        boolean makeAllMove = params.makeAllMove;
        boolean doGraphics = params.doGraphics;
        boolean doElectrostatics = params.doElectrostatics;
        boolean ifMOPPresent = params.ifMOPPresent;
        boolean ifGrapheneDataFile = params.ifGrapheneDataFile;
        boolean ifGasOPLS = params.ifGasOPLS;
        boolean ifCOFPresent = params.ifCOF;
        // double mu1 = -params.mu1*BOLTZMANN_K*Kelvin.UNIT.toSim(temperature);
        boolean isGasCOMPASS = params.isGasCOMPASS;
        boolean isGasTraPPE = params.isGasTraPPE;
        double mu1 = params.mu1;
        double muDecrease = params.muDecrease;
        double multiplier = params.multiplier;
        //mu = -KbTexp(-BU)
        double mu2 = params.mu2;
        double sigma = params.sigma;
        List<Integer> muValues = new ArrayList<>();
        List<String> mu2Values = new ArrayList<>();
        mu2Values.add(params.geomOne);
        List<Integer> temperatureList = new ArrayList<>();

        List<String> mopNames = new ArrayList<>();
       /* mopNames.add(params.confNameGraphene);
        mopNames.add(params.confNameGraphene1);
        mopNames.add(params.confNameGraphene2);
        mopNames.add(params.confNameGraphene3);
        mopNames.add(params.confNameGraphene4);
        mopNames.add(params.confNameGraphene5);
        mopNames.add(params.confNameGraphene6);
        mopNames.add(params.confNameGraphene7);
        mopNames.add(params.confNameGraphene8);
        mopNames.add(params.confNameGraphene9);*/
        mopNames.add(params.confNameGraph1);
       /* mopNames.add(params.confNameGraph2);
        mopNames.add(params.confNameGraph3);
        mopNames.add(params.confNameGraph4);
        mopNames.add(params.confNameGraph5);
        mopNames.add(params.confNameGraph6);
        mopNames.add(params.confNameGraph7);
        mopNames.add(params.confNameGraph8);
        mopNames.add(params.confNameGraph9);
        mopNames.add(params.confNameGraph10);
        mopNames.add(params.confNameGraph11);
        mopNames.add(params.confNameGraph12);*/

      /*  List<String> mopTypes = new ArrayList<>();

        mopTypes.add(params.confNameEight4);*/

        List<double[]> boxLen = new ArrayList<>();
        double[] boxOne = new double[3];
        List<double[]> grapheneSheetsPosnt = new ArrayList<>();
        double[] graphPosn = new double[2];
       /* for (int delta = 0; delta < 1; delta++){
            double val = 5 + delta * 0.5;
            boxOne = new double[]{30, 30, 2*val};
            boxLen.add(boxOne);
            graphPosn = new double[]{-val/2, val/2};
            grapheneSheetsPosnt.add(graphPosn);
        }*/
        boxOne = new double[]{18.5,18.5,50};
        boxLen.add(boxOne);
        System.out.println(Arrays.deepToString(boxLen.toArray()));
        System.out.println(Arrays.deepToString(grapheneSheetsPosnt.toArray()));
        List<String> gasSim = new ArrayList<>();
        gasSim.add(params.confNamegas);
        //gasSim.add(params.confNamegasOne);

        double t1Start = System.nanoTime();
        for (int i = (int) mu1; i>params.muLimit; i= (int) (i+muDecrease)){
            muValues.add(i);
        }

        System.out.println(muValues);
        int num = 0;
        String confNameMOPOne = new String("");

        System.out.println(temperature);

        //String confNameMOPOne = params.confNameEight1;
    //    for(int gamma = 0; gamma< mu2Values.size(); gamma++) {
          /*  String geom = mopTypes.get(gamma);
            System.out.println("geomNames " +geom);*/
            // mu2 = mu2Values.get(gamma);
           // String geom =mu2Values.get(gamma);
        String geom = "tetra";
        for (int phi = 0; phi < mopNames.size(); phi++){
            confNameMOPOne =mopNames.get(phi);
            System.out.println(confNameMOPOne);
            for (int a = 0; a < gasSim.size(); a++) {
                String gasSimAct = gasSim.get(a);
                System.out.println(gasSimAct);
                //String fileOutPutName = filNames.get(a);
                // double sideLength = params.boxSize.getX(0);
                for (int p = 0; p < boxLen.size(); p++) {
                    // double temperature = mopNames.get(p);

                   // double[] levels = grapheneSheetsPosnt.get(p);
                    double[] boxSizeN = boxLen.get(p);
                    double[] intn = boxLen.get(0);
                    Vector boxSizeAct = Vector.of(boxSizeN[0], boxSizeN[1], boxSizeN[2]) ;
                    truncatedRadiusLJ = 12.5;
                    System.out.println("Box " + boxSizeAct+ " "+ truncatedRadiusLJ+ " ");
                    for (int i = 0; i < muValues.size(); i++) {
                        long numSteps = params.numSteps;
                        int mu = muValues.get(i);
                       // System.out.println("mu1 " + mu);
                        //  int mu = -3500;
                        //System.out.println("temperature " + params.temperature);
                        GCMCGO sim = new GCMCGO( gasSimAct, confNameGraphene,  grapheneOne, grapheneTwo, grapheneThree, grapheneFour, numAtomOne,temperature, truncatedRadius, truncatedRadiusLJ, sigma, mu, mu2, ifGraphenePresent, ifMultipleGraphenePresent, ifMoveRotateMoves, boxSizeAct, makeAllMove, doElectrostatics, ifGrapheneDataFile, isGasCOMPASS, isGasTraPPE, ifGasOPLS) ;
                        if (massMembrane < 1) {
                            massMembrane = 1;
                        }

                        double massMOP, massGraphene, massGas;
                        if (ifMOPPresent) {
                            massMOP = sim.speciesMOP.getMass() * numAtomTwo;
                        } else {
                            massMOP = 0;
                        }
                        if (ifGraphenePresent) {
                            massGraphene = sim.speciesGrapheneOne.getMass() * numAtomOne;
                        } else {
                            massGraphene = 0;
                        }

                        if (!ifMOPPresent && !ifGraphenePresent) {
                            massGraphene = 1;
                            massMOP = 1;
                        }
                        //  System.out.println(sim.speciesMOP.getMass());
                        //  System.out.println(confNameGraphene);
                        if (doGraphics) {
                            final SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, "GCMC");
                            DiameterHashByType dhbt = (DiameterHashByType) simGraphic.getDisplayBox(sim.box).getDiameterHash();
                            if (ifMOPPresent) {
                                //  dhbt.setDiameter(sim.speciesMOP.getAtomType(0), 1);
                                // dhbt.setDiameter(sim.speciesMOP.getAtomType(0), 1);
                                // dhbt.setDiameter(sim.speciesMOP.getAtomType(1), 0.8);
                                // dhbt.setDiameter(sim.speciesMOP.getTypeByName("Cu"), 1);
                        /* dhbt.setDiameter(sim.speciesMOP.getAtomType(2), 0.8);
                         dhbt.setDiameter(sim.speciesMOP.getAtomType(3), 0.8);
                         dhbt.setDiameter(sim.speciesMOP.getAtomType(4), 0.8);*/
                                // ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesMOP.getAtomTypes().get(0), Color.green);
                                // ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesMOP.getAtomTypes().get(1), Color.red);
                                // ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesMOP.getAtomTypes().get(2), ColorExtra.goldenRod);
                                // ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesMOP.getAtomTypes().get(3), ColorExtra.gray);
                                // ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesMOP.getAtomTypes().get(4), ColorExtra.red);
                                //  dhbt.setDiameter(sim.speciesMOP.getAtomType(5), 0.8);
                            }

                            if (ifGraphenePresent) {
                               /* for (int j = 0; j < sim.speciesGrapheneOne.getUniqueAtomTypes().size(); i++) {
                                    dhbt.setDiameter(sim.speciesGrapheneOne.getAtomType(j), 1.5);
                                }*/
                                dhbt.setDiameter(sim.speciesGas.getAtomType(0), 1.5);
                             //   dhbt.setDiameter(sim.speciesGrapheneOne.getAtomType(1), 1.5);
                               // dhbt.setDiameter(sim.speciesGrapheneOne.getAtomType(2), 1.5);
                               // dhbt.setDiameter(sim.speciesGrapheneOne.getAtomType(3), 1.5);
                                //dhbt.setDiameter(sim.speciesGas.getAtomType(1), 1);
                                //dhbt.setDiameter(sim.speciesGrapheneOne.getAtomType(0), 1.2);
                                //   dhbt.setDiameter(sim.speciesGrapheneOne.getAtomType(0), 1.2);
                                //  dhbt.setDiameter(sim.speciesGrapheneOne.getAtomType(2), 1.2);
                                //  dhbt.setDiameter(sim.speciesGrapheneOne.getAtomType(3), 0.8);
                            /*    dhbt.setDiameter(sim.speciesGrapheneOne.getAtomType(4), 1.2);
                                dhbt.setDiameter(sim.speciesGrapheneOne.getAtomType(5), 1.2);
                                dhbt.setDiameter(sim.speciesGrapheneOne.getAtomType(6), 1.2);
                                dhbt.setDiameter(sim.speciesGrapheneOne.getAtomType(7), 1.2);
                                dhbt.setDiameter(sim.speciesGrapheneOne.getAtomType(8), 1.2);*/
                               /* ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGrapheneOne.getTypeByName("CX"), Color.darkGray);
                                ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGrapheneOne.getTypeByName("CY"), Color.darkGray);
                                ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGrapheneOne.getTypeByName("CZ"), Color.darkGray);
                                ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGrapheneOne.getTypeByName("C4"), ColorExtra.darkGray);
                                ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGrapheneOne.getTypeByName("OL"), ColorExtra.oliveDrab);
                                ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGrapheneOne.getTypeByName("OJ"), ColorExtra.oliveDrab);
                                ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGrapheneOne.getTypeByName("OK"), ColorExtra.oliveDrab);
                                ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGrapheneOne.getTypeByName("OE"), ColorExtra.oliveDrab);
                                ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGrapheneOne.getTypeByName("HK"), ColorExtra.lightcyan1);
                                ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGrapheneOne.getTypeByName("C4"), Color.darkGray);
                                ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGrapheneOne.getTypeByName("HK"), ColorExtra.Violet);
                                ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGrapheneOne.getTypeByName("OJ"), ColorExtra.gold);
                                ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGrapheneOne.getTypeByName("OK"), ColorExtra.cornflowerblue);*/

                            }

                            if (ifGraphenePresent) {
                                ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGrapheneOne.getAtomType(0), Color.darkGray);
                                /*((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGrapheneOne.getAtomType(1), Color.darkGray);
                                ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGrapheneOne.getAtomType(2), ColorExtra.firebrick);
                                ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGrapheneOne.getAtomType(3), ColorExtra.lightGray);*/
                                //  ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGrapheneOne.getAtomType(1), Color.darkGray);
                                //  ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGrapheneOne.getAtomType(2), ColorExtra.darkRed);
                                // ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGrapheneOne.getAtomType(3), Color.lightGray);
                             /*   ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGrapheneOne.getAtomType(4), Color.gray);
                                ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGrapheneOne.getAtomType(5), ColorExtra.lightcyan);
                                ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGrapheneOne.getAtomType(6), ColorExtra.lightcyan);
                                ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGrapheneOne.getAtomType(7), ColorExtra.lightcyan);*/
                            }


                            //MOP
                            if (ifMOPPresent) {
                                //  dhbt.setDiameter(sim.speciesMOP.getAtomType(0), 3);
                  /*  dhbt.setDiameter(sim.speciesMOP.getAtomType(1), 2);
                    dhbt.setDiameter(sim.speciesMOP.getAtomType(2), 1);
                    dhbt.setDiameter(sim.speciesMOP.getAtomType(3), 2);
                    dhbt.setDiameter(sim.speciesMOP.getAtomType(4), 2);
                    dhbt.setDiameter(sim.speciesGas.getAtomType(0), 1.5);*/
                                // dhbt.setDiameter(sim.speciesGas.getAtomType(1), 1);
                                // ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesMOP.getAtomType(0), ColorExtra.copper);
                                //((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesMOP.getAtomType(1), ColorExtra.red);
                                // ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesMOP.getAtomType(2), ColorExtra.red);
                                // ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesMOP.getAtomType(3), ColorExtra.red);
                                //   ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesMOP.getAtomType(4), ColorExtra.lightcyan);
                                // ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesMOP.getAtomType(5), ColorExtra.darkGray);
                                //   ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesMOP.getAtomType(6), ColorExtra.red);
                                //   ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesMOP.getAtomType(7), ColorExtra.blue1);
                                //   ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesMOP.getAtomType(8), ColorExtra.lightcyan);
                                //   ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGas.getAtomType(0), ColorExtra.tomato);
                                //((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGas.getAtomType(1), ColorExtra.lightseagreen);
                            }

                            if (ifSecondGasPresent) dhbt.setDiameter(sim.speciesSecondGas.getAtomType(0), 1.5);
                            //if(ifSecondGasPresent)dhbt.setDiameter(sim.speciesSecondGas.getAtomType(1), 1);

                            //    dhbt.setDiameter(sim.speciesGas.getAtomType(0), 3.5);
                            // dhbt.setDiameter(sim.speciesGas.getAtomType(1), 3);
                            ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGas.getAtomType(0), ColorExtra.lightseagreen);
                            // ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGas.getAtomType(1), ColorExtra.salmon);
                            simGraphic.makeAndDisplayFrame("GCMC");
                            ActivityIntegrate ai2 = new ActivityIntegrate(sim.integrator);
                            sim.getController().addActivity(ai2, Long.MAX_VALUE, 1.0);

                            return;
                        }
                        //System.out.println(sim.box.getMoleculeList(sim.speciesGrapheneOne).size()* sim.speciesGrapheneOne.getMass()+ " " + sim.box.getMoleculeList(sim.speciesGrapheneOne).size()+ " " + sim.speciesGrapheneOne.getMass());

                        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, params.numSteps / 10));
                        int pInterval = 2;
                        long bs = params.numSteps / (pInterval * 5);
                        if (bs == 0) bs = 1;
                        /*MeterPressure pMeter = new MeterPressure(sim.box, sim.integrator.getPotentialCompute());
                        pMeter.setTemperature(sim.integrator.getTemperature());
                        AccumulatorAverage pAccumulator = new AccumulatorAverageFixed(bs);
                        DataPumpListener pPump = new DataPumpListener(pMeter, pAccumulator, pInterval);
                        sim.integrator.getEventManager().addListener(pPump);*/

                        bs = params.numSteps / 50;
                        if (bs == 0) bs = 1;
                        MeterPotentialEnergyFromIntegrator energyMeter = new MeterPotentialEnergyFromIntegrator(sim.integrator);
                        AccumulatorAverage energyAccumulator = new AccumulatorAverageFixed(bs);
                        DataPumpListener uPump = new DataPumpListener(energyMeter, energyAccumulator);
                        sim.integrator.getEventManager().addListener(uPump);

                        nMoleculesGasOne = new MeterNMolecules();
                        nMoleculesGasOne.setBox(sim.box);
                        nMoleculesGasOne.setSpecies(sim.speciesGas);

                        if (ifSecondGasPresent) {
                            nMoleculesGasTwo = new MeterNMolecules();
                            nMoleculesGasTwo.setBox(sim.box);
                            nMoleculesGasTwo.setSpecies(sim.speciesSecondGas);
                        }
                        MeterDensity densityMeter = new MeterDensity(sim.box);
                        AccumulatorAverage densityAccumulator = new AccumulatorAverageFixed(bs);
                        densityMeter.setSpecies(sim.speciesGas);
                        DataPumpListener pumpDensity = new DataPumpListener(densityMeter, densityAccumulator);
                        sim.integrator.getEventManager().addListener(pumpDensity);
                        MeterDensity densityMeterGasTwo = new MeterDensity(sim.box);
                        AccumulatorAverage densityAccumulatorGastwo = new AccumulatorAverageFixed(bs);
                        if (ifSecondGasPresent) {
                            densityMeterGasTwo.setSpecies(sim.speciesSecondGas);
                            DataPumpListener pumpDensityGasTwo = new DataPumpListener(densityMeterGasTwo, densityAccumulatorGastwo);
                            sim.integrator.getEventManager().addListener(pumpDensityGasTwo);
                        }
                        long t1 = System.currentTimeMillis();
                        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, params.numSteps));
                        long t2 = System.currentTimeMillis();
                        //System.out.println("mu " + mu);
                        //  System.out.println("runtime: " + (t2 - t1) * 0.001);
                        Unit pUnit = Bar.UNIT;
                        //     Unit MegaPascal = new PrefixedUnit(Prefix.MEGA, Pascal.UNIT);
                     /*   Unit kiloPascal = new PrefixedUnit(Prefix.KILO, Pascal.UNIT);
                        double avgP = pAccumulator.getData(pAccumulator.AVERAGE).getValue(0);
                        double errP = pAccumulator.getData(pAccumulator.ERROR).getValue(0);
                        double corP = pAccumulator.getData(pAccumulator.BLOCK_CORRELATION).getValue(0);*/
                       // System.out.println("P (Bar) " + pUnit.fromSim(avgP) + " bar " + pUnit.fromSim(errP) + " " + corP);
                        //   System.out.println("P (MPa) " + MegaPascal.fromSim(avgP) + " MPa " + MegaPascal.fromSim(errP) + " " + corP);
                        //  System.out.println("P (kPa) " + kiloPascal.fromSim(avgP) + " kPa " + kiloPascal.fromSim(errP) + " " + corP);
                        //double muR = Constants.BOLTZMANN_K * temperature * Math.log(avgP/temperature);
                        //System.out.println("muR: " + muR);*/

                        double avgRho = densityAccumulator.getData(densityAccumulator.AVERAGE).getValue(0);
                        double errRho = densityAccumulator.getData(densityAccumulator.ERROR).getValue(0);
                        double corRho = densityAccumulator.getData(densityAccumulator.BLOCK_CORRELATION).getValue(0);
                        Unit dm = new PrefixedUnit(Prefix.DECI, Meter.UNIT);
                        Unit dm3 = new CompoundUnit(new Unit[]{dm}, new double[]{3});
                        Unit m3 = new CompoundUnit(new Unit[]{Meter.UNIT}, new double[]{3});
                        Unit kgm3 = new UnitRatio(new PrefixedUnit(Prefix.KILO, Gram.UNIT), m3);
                        Unit kJmol = new UnitRatio(new PrefixedUnit(Prefix.KILO, Joule.UNIT), Mole.UNIT);
                        //Unit moldm3 = new UnitRatio(Mole.UNIT, Liter.UNIT);
                        Unit cm = new PrefixedUnit(Prefix.CENTI, Meter.UNIT);
                        Unit cm3 = new CompoundUnit(new Unit[]{cm}, new double[]{3});
                        Unit nm = new PrefixedUnit(Prefix.NANO, Meter.UNIT);
                        Unit nm3 = new CompoundUnit(new Unit[]{nm}, new double[]{3});
                        //  Unit gcm3 = new UnitRatio(Gram.UNIT, cm3);
                        //   Unit cm3g = new UnitRatio(cm3, Gram.UNIT);
                        double numAtomsAvg = avgRho * sim.box.getBoundary().volume();
                        double errnum = errRho * sim.box.getBoundary().volume();
                        //    System.out.println(sim.box.getBoundary().volume()+ " volume");
                        double cornum = corRho;
                        Unit kilogram = new PrefixedUnit(Prefix.KILO, Gram.UNIT);
                        // double mopMass = sim.speciesMOP.getMass();
                        //  System.out.println(sim.speciesGrapheneOne.getMass());
                        double grapheneMass = sim.speciesGrapheneOne.getMass();
                        double mmolPg = (numAtomsAvg* 1000)/(grapheneMass * 4);
                        //double errmmolPg = (errnum* 1000)/mopMass;
                        //Unit moleperKilo = new UnitRatio(Mole.UNIT, kilogram);
                        double err1 = errnum *100 / numAtomsAvg;
                        System.out.println("Num atoms One : " + numAtomsAvg + " " + errnum  +" "+ err1  +" "+ mmolPg);
                        // double massPer = (numAtomsAvg* sim.speciesGas.getMass())/((sim.speciesMOP.getMass())+ (numAtomsAvg* sim.speciesGas.getMass()));
                        // System.out.println("mass Per " +massPer*100);
                        //System.out.println("mol/mass : "+ numAtomsAvg/ sim.speciesMOP.getMass() + " "+ sim.speciesMOP.getMass() + " " + moleperKilo.fromSim(numAtomsAvg/ sim.speciesMOP.getMass()));
                        //  System.out.println("Actual "+ avgRho);
                        //System.out.println("nm " + nm3.fromSim(sim.box.getBoundary().volume()));
                        //System.out.println("Num atoms : " + numAtomsAvg + " " + errnum + " " + cornum + " " + (sim.speciesGas.getMass()*numAtomsAvg*1000/sim.speciesMOP.getMass()));
                        // System.out.println("num atoms (/cm3) "+cm3.fromSim( sim.box.getBoundary().volume())*errRho + " " + cm3.fromSim( sim.box.getBoundary().volume())*errRho);
                        // System.out.println("num atoms (/m3) " +m3.fromSim( sim.box.getBoundary().volume())*errRho + " " + m3.fromSim( sim.box.getBoundary().volume())*errRho);
                        // System.out.println("num atoms (/nm3) " +nm3.fromSim( sim.box.getBoundary().volume())*errRho + " " + nm3.fromSim( sim.box.getBoundary().volume())*errRho);
                        // System.out.println("num atoms (/nm3) " +nm_3.fromSim(avgRho) + " " + nm_3.fromSim(errRho) );
                        // System.out.println("Rho (mol/dm3) "+kgm3.fromSim(avgRho)/sim.speciesGas.getMass() + " " + sim.speciesGas.getMass());
                        //AtomType C1 = new AtomType(Carbon.INSTANCE, "C1");
                        //System.out.println("Amount sorbed " + numAtomsAvg*1000*sim.speciesGas.getMass() /(sim.speciesMOP.getMass()) + " err " +errnum*1000*sim.speciesGas.getMass() /(sim.speciesMOP.getMass()));


                        // massGas = numAtomsAvg * sim.speciesGas.getMass();
                        // System.out.println(massGas + " Graphene " + massGraphene  + " single graphene " + sim.speciesGrapheneOne.getMass());
                        //  System.out.println("Mass ratio " + massGas/(massGraphene));
                        //   System.out.println(massMOP + " "+massGraphene );
                        //  System.out.println("Amount sorbed " + 1000*massGas/(massGraphene+massMOP)  + " err " + 1000*errnum*sim.speciesGas.getMass()/(massGraphene+massMOP));
                        //      System.out.println("Amount sorbed per Graphene basis " + massGas/((massGraphene))  + " err " +errnum*sim.speciesGas.getMass()/(massGraphene));
                        //  System.out.println("Density (g/cm3) " + gcm3.fromSim(avgRho));
                        //System.out.println("Mass fraction : "+sim.speciesGas.getMass()*1000/(sim.speciesMOP.getMass()));
          /*  double avgPE = energyAccumulator.getData(energyAccumulator.AVERAGE).getValue(0) / numAtomsAvg;
            double errPE = energyAccumulator.getData(energyAccumulator.ERROR).getValue(0) / numAtomsAvg;
            double corPE = energyAccumulator.getData(energyAccumulator.BLOCK_CORRELATION).getValue(0);*/
                        //  System.out.println("PE " + avgPE + " " + errPE + " " + corPE);
                        // double numAtomGasOne = nMoleculesGasOne.getDataAsScalar();
                        //System.out.println("Atoms gas One: " +numAtomGasOne);

                        if (ifSecondGasPresent) {
                            double avgRho2 = densityAccumulatorGastwo.getData(densityAccumulator.AVERAGE).getValue(0);
                            double errRho2 = densityAccumulatorGastwo.getData(densityAccumulator.ERROR).getValue(0);
                            double corRho2 = densityAccumulatorGastwo.getData(densityAccumulator.BLOCK_CORRELATION).getValue(0);
                            double numAtomsAvg2 = avgRho2 * sim.box.getBoundary().volume();
                            double errnum2 = errRho2 * sim.box.getBoundary().volume();
                            double err2 = errnum2 * 100 / numAtomsAvg2;
                            //System.out.println("Num atoms : " + numAtomsAvg +" "+numAtomsAvg2);
                            //  System.out.println("Error " + err1 +" "+ err2+" \n");
                        }

//        System.out.println(numAtomsAvg*sim.speciesGas.getMass()*100/(sim.speciesGrapheneOne.getMass()) + " " +numAtomsAvg*sim.speciesGas.getMass()+ " "+sim.speciesGrapheneOne.getMass());
                        //System.out.println("rho (kg/m3) " + kgm3.fromSim(avgRho) + " " + kgm3.fromSim(errRho) + " " + corRho);
                        //  System.out.println("rho " + (avgRho) + " " + errRho + " " + corRho + " " + (errRho / avgRho) * 100 );
                        //System.out.println("\n");
                        // System.out.println("uptake  " + 22.414 * kgm3.fromSim(avgRho) / sim.speciesGas.getMass());
                        // System.out.println("Rho (mol/dm3) " + moldm3.fromSim(avgRho) + " " +errRho + " " +corRho);
                        // expected values based on 10^8 steps
                        // stdev based on 10^8 steps for N=500 uncertainty, scaled up to 10^6 steps
                        //   other sims have significant correlation
                        // 4 sigma should fail 1 in 16,000 runs;
                    }
                }
                double t2 = System.nanoTime();
                System.out.println((t2 - t1Start) / Math.pow(10, 9));
            }
        }
    }

    public Integrator getIntegrator(){
        return integrator;
    }
    public static class GCMCGOParams extends ParameterBase {
        public Vector grapheneThirteen = new Vector3D(15.0,15.0, 60.0);
        public Vector grapheneOne = new Vector3D(0.0,0.0, 4.0);
        public Vector grapheneTwo = new Vector3D(0,0, -4.0);
        public Vector grapheneThree = new Vector3D(0,0, 12.0);
        public Vector grapheneFour = new Vector3D(0,0, -12.0);
        public Vector grapheneSeven = new Vector3D(0,0, 38.5);
        public Vector grapheneSix = new Vector3D(0,0, -19.5);

        public Vector grapheneNine = new Vector3D(15.0,15.0, 40.0);
        public Vector grapheneTen = new Vector3D(15.0,15.0, 50.0);
        public Vector grapheneFive = new Vector3D(15.0,15.0, -10.0);
        public Vector grapheneEight = new Vector3D(15.0,15.0, -40.0);
        public Vector grapheneEleven = new Vector3D(15.0,15.0, -50.0);
        public Vector grapheneTwelve = new Vector3D(15.0,15.0, -60.0);
        public Vector centreMOPTwo = new Vector3D(0.0,0.0,30);
        public Vector centreMOP = new Vector3D(0.0,0.0,-30);
        public int numAtomOne = 1;
        public int numAtomTwo = 1;

        public double sigma =2.7;
        public boolean makeAllMove = false;
        public boolean ifCOF = false;
        public String confNameCOF = "F://Avagadro//mop//tetra_CuCOOCH3";
        public double boxZLen = 30;
        public double graphenezOne = -6.25;
        public double graphenezTwo = 6.25;
        public double graphenezThree = 18.75;
        public double graphenezFour = -18.75;
        public String confNameGraphene = "D:\\Sem-X\\GO\\graphitis\\GO_sheet";
        public boolean ifGasfileRead = false;
        public String confNameGraphene1 = "D:\\Sem-IX\\GO\\go0019";
        public String confNameGraphene2 = "D:\\Sem-IX\\GO\\go0038";
        public String confNameGraphene3 = "D:\\Sem-IX\\GO\\go0056";
        public String confNameGraphene4 = "D:\\Sem-IX\\GO\\go0075";
        public String confNameGraphene5 = "D:\\Sem-IX\\GO\\go0094";
        public String confNameGraphene6 = "D:\\Sem-IX\\GO\\go00113";
        public String confNameGraphene7 = "D:\\Sem-IX\\GO\\go00132";
        public String confNameGraphene8 = "D:\\Sem-IX\\GO\\go00150";
        public String confNameGraphene9 = "D:\\Sem-IX\\GO\\go00169";
        public String confNameGraph1 ="D:\\Sem-IX\\GO\\allSheets\\acid\\go10COOH";
        public String confNameGraph2 ="D:\\Sem-IX\\GO\\allSheets\\acid\\go20COOH";
        public String confNameGraph3 ="D:\\Sem-IX\\GO\\allSheets\\acid\\go29COOH";
        public String confNameGraph4 ="D:\\Sem-IX\\GO\\allSheets\\acid\\go38COOH";
        public String confNameGraph5 ="D:\\Sem-IX\\GO\\allSheets\\acid\\go47COOH";
        public String confNameGraph6 ="D:\\Sem-IX\\GO\\allSheets\\ether\\go019eth";
        public String confNameGraph7 ="D:\\Sem-IX\\GO\\allSheets\\ether\\go038eth";
        public String confNameGraph8 ="D:\\Sem-IX\\GO\\allSheets\\ether\\go056eth";
        public String confNameGraph9 ="D:\\Sem-IX\\GO\\allSheets\\ether\\go075eth";
        public String confNameGraph10 ="D:\\Sem-IX\\GO\\allSheets\\ether\\go094eth";
        public String confNameGraph11 ="D:\\Sem-IX\\GO\\allSheets\\ether\\go0113eth";
        public String confNameGraph12 ="D:\\Sem-IX\\GO\\allSheets\\ether\\go0132eth";
        public Vector boxSize = new Vector3D(16, 16, 16);
        public double temperature = 298;
        public double mu2 = -470;
        public double mu1 = -1400;
        public double muDecrease = -20;
        public int muLimit = -3200;
        public long numSteps =10000;
        // public double side =50.0;
        public boolean ifXYZfile = false;

        public String confNamegas = "D:\\Sem-X\\GO\\gasOPLS\\ethane";
        public String confNamegasOne = "ch4" ;
        public String confNameGasTwo = "ch4" ;
        public String confNameGasThree = "methylpropane" ;
        public String confNamegasFour = "methylpropene" ;
        public String confNamegasFive = "cisbutene" ;
        public String confNamegasSix = "transbutene" ;
        public String confNamegasSeven = "ethane" ;
        public String confNamegasEight = "ethene" ;
        public String confNamegasNine = "ch4" ;
        public String confNamegasTen = "propane" ;
        public String confNamegasEleven = "propene" ;
        public boolean isGasTraPPE = false;
        public boolean isGasOPLS = true;
        public String geomOne = "icosa" ;
        public double truncatedRadiusLJ =12.8;
        public double truncatedRadius = truncatedRadiusLJ;
        public String struc = "edge";
        public boolean doFFAnalaysis = true;
        public boolean doGraphics = false;
        public boolean ifautoMOP = false;
        public boolean ifGasOPLS = true;
        public boolean isGasCOMPASS = false;
        public boolean ifGraphenePresent = false;
        public boolean ifGrapheneDataFile = true;
        public boolean ifMultipleGraphenePresent = true;
        public boolean ifSecondGasPresent = false;
        public boolean ifMOPPresent = false;
        public boolean ifMoveRotateMoves = true;
        public boolean doElectrostatics = false;
        public boolean doMOPPDBOutput = false;
        public double multiplier = 9;
    }
    public List<List<AtomType>> pairsAtomsSecondGas = new ArrayList<>();
    public ArrayList<ArrayList<Integer>> connectedAtomsSecondGas, connectivityModifiedSecondGas = new ArrayList<>();
    public Map<Integer,String> atomMapSecondGas = new HashMap<>();
    public ArrayList<Integer> bondListSecondGas, bondsNumSecondGas = new ArrayList<>();
    public Map<String, double[]> atomicPotMapSecondGas = new HashMap<>();
    public Map<Integer, String> atomIdentifierMapModifiedSecondGas = new HashMap<>();
    public List<int[]> dupletsSortedSecondGas, tripletsSortedSecondGas, quadrupletsSortedSecondGas = new ArrayList<>();
    public Map<String[],List<int[]>> bondTypesMapSecondGas, angleTypesMapSecondGas, torsionTypesMapSecondGas = new HashMap<>();
    public static MeterNMolecules nMoleculesGasTwo, nMoleculesGasOne;
}

