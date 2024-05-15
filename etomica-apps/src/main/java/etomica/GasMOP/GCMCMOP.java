package etomica.GasMOP;

import etomica.action.BoxInflate;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.DiameterHashByType;
import etomica.box.Box;
import etomica.chem.elements.Carbon;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.DataPumpListener;
import etomica.data.meter.MeterDensity;
import etomica.data.meter.MeterNMolecules;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.meter.MeterPressure;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.Integrator;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.*;
import etomica.molecule.IMolecule;
import etomica.molecule.MoleculeSourceRandomMolecule;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.potential.*;
import etomica.potential.COMPASS.LJCOMPASS;
import etomica.potential.COMPASS.PDBReaderCOMPASS;
import etomica.potential.TraPPE.speciesGasTraPPE;
import etomica.potential.UFF.*;
import etomica.simulation.Simulation;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesManager;
import etomica.units.*;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.util.random.RandomMersenneTwister;
import org.apache.bcel.Constants;

import java.awt.*;
import java.util.*;
import java.util.List;

import static etomica.potential.TraPPE.speciesGasTraPPE.ChemForm.*;
import static etomica.util.Constants.BOLTZMANN_K;

public class GCMCMOP extends Simulation {
    public IntegratorMC integrator;
    public MCMoveMolecule mcMoveMolecule;
    public MCMoveMoleculeRotate mcMoveMoleculeRotate;
    public MCMoveInsertDelete mcMoveID, mcMoveIDTwo;
    public ISpecies speciesMOP, speciesGas, speciesGrapheneOne, speciesGrapheneTwo, speciesGrapheneThree, speciesGrapheneFour,speciesGrapheneFive,speciesGrapheneSix,speciesGrapheneSeven,speciesGrapheneEight, speciesSecondGas,speciesGrapheneNine,speciesGrapheneTen,speciesGrapheneEleven,speciesGrapheneTwelve, speciesGrapheneThirteen;
    public ISpecies speciesMOPTwo, speciesMOPThree,speciesMOPFour, speciesMOPFive;
    public Box box;
    public P2LennardJones potential;
    public SpeciesManager sm;
    public static int atom1, atom2, atom3;
    public int numCarbon;
    public static double sigmaIKey, sigmaJKey, epsilonIKey, epsilonJKey;
    public int i;
    public static String atomName1,atomName2, atomName3 ;
    public static List<List<AtomType>> listMOPGasMixed = new ArrayList<>();
    public static List<List<AtomType>> listGrapheneGasMixed = new ArrayList<>();
    public GCMCMOP(String confNameOne, String confNameGasOne, String confNameGasTwo, String confNameGraphene, Vector centreMOP, Vector centreMOPTwo, Vector centreMOPThree,  Vector centreMOPFour, Vector grapheneOne, Vector grapheneTwo,Vector grapheneThree,Vector grapheneFour,Vector grapheneFive,Vector grapheneSix,Vector grapheneSeven,Vector grapheneEight,Vector grapheneNine, Vector grapheneTen, Vector grapheneEleven, Vector grapheneTwelve, Vector grapheneThirteen, int numMolOne, double temperature, double truncatedRadius, double truncatedRadiusLJ, double sigma, double mu1, double mu2, boolean ifGraphenePresent, boolean ifSecondGasPresent, boolean ifMultipleGraphenePresent, boolean ifMoveRotateMoves, Vector boxSize, boolean makeAllMove, boolean doElectrostatics, double multiplier, boolean isGasCOMPASS, boolean isGasTraPPE){
        super(Space3D.getInstance());

        //Make Species
       // setRandom(new RandomMersenneTwister(random));
        PDBReaderMOP pdbReaderMOP = new PDBReaderMOP();
        PDBReaderMOP pdbReaderMOP2 = new PDBReaderMOP();
        PDBReaderReplica pdbReaderReplica = new PDBReaderReplica();
        PDBReaderCOMPASS pdbReaderCOMPASS = new PDBReaderCOMPASS();
        PDBReaderReplicaNew pdbReaderReplicaNew = new PDBReaderReplicaNew();
        GeneralGrapheneReader grapheneReader = new GeneralGrapheneReader();
        SetPotential SetPotential = new SetPotential();


        if(ifGraphenePresent){
            if(ifMultipleGraphenePresent){
                speciesGrapheneOne = grapheneReader.getSpecies(confNameGraphene, grapheneOne, false);
                addSpecies(speciesGrapheneOne);
               /* speciesGrapheneTwo = grapheneReader.getSpecies(confNameGraphene, grapheneTwo, false);
               // speciesGrapheneThree = grapheneReader.getSpecies(confNameGraphene, grapheneThree, false);
                speciesGrapheneFour = grapheneReader.getSpecies(confNameGraphene, grapheneFour, false);
                speciesGrapheneFive = grapheneReader.getSpecies(confNameGraphene, grapheneFive, false);
                //speciesGrapheneSix = grapheneReader.getSpecies(confNameGraphene, grapheneSix, false);
                speciesGrapheneSeven = grapheneReader.getSpecies(confNameGraphene, grapheneSeven, false);
                speciesGrapheneEight = grapheneReader.getSpecies(confNameGraphene, grapheneEight, false);
                speciesGrapheneNine = grapheneReader.getSpecies(confNameGraphene, grapheneNine, false);
             //   speciesGrapheneTen = grapheneReader.getSpecies(confNameGraphene, grapheneTen, false);
              //  speciesGrapheneEleven = grapheneReader.getSpecies(confNameGraphene, grapheneEleven, false);
                speciesGrapheneTwelve = grapheneReader.getSpecies(confNameGraphene, grapheneTwelve, false);
                speciesGrapheneThirteen = grapheneReader.getSpecies(confNameGraphene, grapheneThirteen, false);
               // addSpecies(speciesGrapheneOne);
                addSpecies(speciesGrapheneTwo);
               // addSpecies(speciesGrapheneThree);
                addSpecies(speciesGrapheneFour);
                addSpecies(speciesGrapheneFive);
               // addSpecies(speciesGrapheneSix);
                addSpecies(speciesGrapheneSeven);
                addSpecies(speciesGrapheneEight);
                addSpecies(speciesGrapheneNine);
                //addSpecies(speciesGrapheneTen);
               // addSpecies(speciesGrapheneEleven);
                addSpecies(speciesGrapheneTwelve);
                addSpecies(speciesGrapheneThirteen);
                numCarbon = grapheneReader.getNumCarbon(confNameGraphene);
                System.out.println(numCarbon);*/
            }else {
                speciesGrapheneOne = grapheneReader.getSpecies(confNameGraphene, grapheneOne, false);
                speciesGrapheneFive = grapheneReader.getSpecies(confNameGraphene, grapheneFive, false);
                addSpecies(speciesGrapheneOne);
                addSpecies(speciesGrapheneFive);
            }
        }
        speciesMOP = pdbReaderMOP.getSpeciesMOP(confNameOne, false, new Vector3D(0.0,0.0,0.0), false);
        System.out.println(speciesMOP.getMass());
        cifReader cifReader = new cifReader();
        speciesGasTraPPE speciesGasTraPPE = new speciesGasTraPPE();
        if(isGasTraPPE){
            if(confNameGasOne.equals("F://Avagadro//molecule//ch4")){
                etomica.potential.TraPPE.speciesGasTraPPE.ChemForm = CH4;
            } else if (confNameGasOne.equals("F://Avagadro//molecule//ethane")) {
                etomica.potential.TraPPE.speciesGasTraPPE.ChemForm = C2H6;
            } else if (confNameGasOne.equals("F://Avagadro//molecule//propane")) {
                etomica.potential.TraPPE.speciesGasTraPPE.ChemForm = C3H8;
            } else if (confNameGasOne.equals("F://Avagadro//molecule//ethene")) {
                etomica.potential.TraPPE.speciesGasTraPPE.ChemForm = C2H4;
            } else if (confNameGasOne.equals("F://Avagadro//molecule//propene")) {
                etomica.potential.TraPPE.speciesGasTraPPE.ChemForm = C3H6;
            }
        }
        //speciesMOP = cifReader.speciesCIF(confNameOne,false);
        addSpecies(speciesMOP);
        //species Gas
        if(isGasCOMPASS){
            speciesGas = pdbReaderCOMPASS.getSpecies(confNameGasOne, false, true, new Vector3D(0,0,0));
        }else if (isGasTraPPE){
            speciesGas = speciesGasTraPPE.speciesGasTraPPE(Space3D.getInstance(), etomica.potential.TraPPE.speciesGasTraPPE.ChemForm, false);
        } else {
            speciesGas = pdbReaderReplica.getSpecies(confNameGasOne, true,new Vector3D(0,0,0), false);
        }
        addSpecies(speciesGas);

        //SecondGas and Bonding Parameters of Second Gas
      if(ifSecondGasPresent){
            speciesSecondGas = pdbReaderReplicaNew.getSpecies(confNameGasTwo);
            addSpecies(speciesSecondGas);
            pairsAtomsSecondGas = SetPotential.getSpeciesPairs(speciesSecondGas);
            connectedAtomsSecondGas =pdbReaderReplicaNew.getConnectivityWithoutRunning();
            connectivityModifiedSecondGas = pdbReaderReplicaNew.getConnectivityModifiedWithoutRunning();
            atomMapSecondGas = pdbReaderReplicaNew.getAtomMapWithoutRunning();
            bondListSecondGas = pdbReaderReplicaNew.getBondList(connectedAtomsSecondGas, atomMapSecondGas);
            atomicPotMapSecondGas = pdbReaderReplicaNew.atomicPotMap();
            atomIdentifierMapModifiedSecondGas = pdbReaderReplicaNew.getatomIdentifierMapModified();
            dupletsSortedSecondGas= pdbReaderReplicaNew.getDupletesSorted();
            tripletsSortedSecondGas=pdbReaderReplicaNew.getAnglesSorted();
            quadrupletsSortedSecondGas=pdbReaderReplicaNew.getTorsionSorted();
            bondsNumSecondGas = pdbReaderReplicaNew.getBonds();
            bondTypesMapSecondGas= pdbReaderReplicaNew.idenBondTypes(dupletsSortedSecondGas, atomIdentifierMapModifiedSecondGas);
            angleTypesMapSecondGas= pdbReaderReplicaNew.idenAngleTypes(tripletsSortedSecondGas, atomIdentifierMapModifiedSecondGas);
            torsionTypesMapSecondGas= pdbReaderReplicaNew.idenTorsionTypes(quadrupletsSortedSecondGas, atomIdentifierMapModifiedSecondGas);
        }

      box = this.makeBox();

      box.getBoundary().setBoxSize(boxSize);
      List<Vector> oldPositions = new ArrayList<>();
    /*  box.setNMolecules(speciesMOP, 5);
      List<Vector> oldPositions = new ArrayList<>();
        IMolecule moleculeMOPOne = box.getMoleculeList().get(0);
      while (oldPositions.size() < moleculeMOPOne.getChildList().size()) {
          oldPositions.add(space.makeVector());
      }
        Vector originTwo = new Vector3D(0,0, 50);
        IMolecule moleculeMOPTwo = box.getMoleculeList().get(1);
      moleculeMOPTwo.getChildList().forEach(atom -> {
          oldPositions.get(atom.getIndex()).E(atom.getPosition());
          atom.getPosition().PE(originTwo);
          Vector shift = box.getBoundary().centralImage(atom.getPosition());
          atom.getPosition().PE(shift);
      });
      Vector originOne = new Vector3D(0,0, 20);
        moleculeMOPOne.getChildList().forEach(atom -> {
            oldPositions.get(atom.getIndex()).E(atom.getPosition());
            atom.getPosition().PE(originOne);
            Vector shift = box.getBoundary().centralImage(atom.getPosition());
            atom.getPosition().PE(shift);
        });
        IMolecule moleculeMOPThree = box.getMoleculeList().get(2);
        Vector originThree = new Vector3D(0,0, 0);
        moleculeMOPThree.getChildList().forEach(atom -> {
            oldPositions.get(atom.getIndex()).E(atom.getPosition());
            atom.getPosition().PE(originThree);
            Vector shift = box.getBoundary().centralImage(atom.getPosition());
            atom.getPosition().PE(shift);
        });
        IMolecule moleculeMOPFour = box.getMoleculeList().get(3);
        Vector originFour = new Vector3D(0,0,-20 );
        moleculeMOPFour.getChildList().forEach(atom -> {
            oldPositions.get(atom.getIndex()).E(atom.getPosition());
            atom.getPosition().PE(originFour);
            Vector shift = box.getBoundary().centralImage(atom.getPosition());
            atom.getPosition().PE(shift);
        });
        IMolecule moleculeMOPFive = box.getMoleculeList().get(4);
        Vector originFive = new Vector3D(0,0,-50 );
        moleculeMOPFive.getChildList().forEach(atom -> {
            oldPositions.get(atom.getIndex()).E(atom.getPosition());
            atom.getPosition().PE(originFive);
            Vector shift = box.getBoundary().centralImage(atom.getPosition());
            atom.getPosition().PE(shift);
        });*/
        if(ifGraphenePresent ){
            if(ifMultipleGraphenePresent ){
                box.setNMolecules(speciesGrapheneOne, 2);
                IMolecule moleculeMOPZero = box.getMoleculeList().get(0);
                while (oldPositions.size() < moleculeMOPZero.getChildList().size()) {
                    oldPositions.add(space.makeVector());
                }
             /*   Vector originZero = new Vector3D(0,0, 50);
                moleculeMOPZero.getChildList().forEach(atom -> {
                    oldPositions.get(atom.getIndex()).E(atom.getPosition());
                    atom.getPosition().PE(originZero);
                    Vector shift = box.getBoundary().centralImage(atom.getPosition());
                    atom.getPosition().PE(shift);
                });
                IMolecule moleculeMOPOne = box.getMoleculeList().get(1);
                Vector originOne = new Vector3D(0,0, 40);
                moleculeMOPOne.getChildList().forEach(atom -> {
                    oldPositions.get(atom.getIndex()).E(atom.getPosition());
                    atom.getPosition().PE(originOne);
                    Vector shift = box.getBoundary().centralImage(atom.getPosition());
                    atom.getPosition().PE(shift);
                });*/
              /*  IMolecule moleculeMOPThree = box.getMoleculeList().get(1);
                Vector originThree = new Vector3D(0,0, 30);
                moleculeMOPThree.getChildList().forEach(atom -> {
                    oldPositions.get(atom.getIndex()).E(atom.getPosition());
                    atom.getPosition().PE(originThree);
                    Vector shift = box.getBoundary().centralImage(atom.getPosition());
                    atom.getPosition().PE(shift);
                });
                IMolecule moleculeMOPFour = box.getMoleculeList().get(2);
                Vector originFour = new Vector3D(0,0, 20);
                moleculeMOPFour.getChildList().forEach(atom -> {
                    oldPositions.get(atom.getIndex()).E(atom.getPosition());
                    atom.getPosition().PE(originFour);
                    Vector shift = box.getBoundary().centralImage(atom.getPosition());
                    atom.getPosition().PE(shift);
                });*/
                IMolecule moleculeMOPFive = box.getMoleculeList().get(0);
                moleculeMOPFive.getChildList().forEach(atom -> {
                    oldPositions.get(atom.getIndex()).E(atom.getPosition());
                    atom.getPosition().PE(grapheneTwo);
                    Vector shift = box.getBoundary().centralImage(atom.getPosition());
                    atom.getPosition().PE(shift);
                });
                IMolecule moleculeMOPSix = box.getMoleculeList().get(1);
                moleculeMOPSix.getChildList().forEach(atom -> {
                    oldPositions.get(atom.getIndex()).E(atom.getPosition());
                    atom.getPosition().PE(grapheneSix);
                    Vector shift = box.getBoundary().centralImage(atom.getPosition());
                    atom.getPosition().PE(shift);
                });
              /*  IMolecule moleculeMOPSeven = box.getMoleculeList().get(5);
                Vector originSeven = new Vector3D(0,0, -20);
                moleculeMOPSeven.getChildList().forEach(atom -> {
                    oldPositions.get(atom.getIndex()).E(atom.getPosition());
                    atom.getPosition().PE(originSeven);
                    Vector shift = box.getBoundary().centralImage(atom.getPosition());
                    atom.getPosition().PE(shift);
                });
                IMolecule moleculeMOPEight = box.getMoleculeList().get(0);
                Vector originEight = new Vector3D(0,0, -30);
                moleculeMOPEight.getChildList().forEach(atom -> {
                    oldPositions.get(atom.getIndex()).E(atom.getPosition());
                    atom.getPosition().PE(originEight);
                    Vector shift = box.getBoundary().centralImage(atom.getPosition());
                    atom.getPosition().PE(shift);
                });*/
             /*   IMolecule moleculeMOPNine = box.getMoleculeList().get(8);
                Vector originNine = new Vector3D(0,0, -40);
                moleculeMOPNine.getChildList().forEach(atom -> {
                    oldPositions.get(atom.getIndex()).E(atom.getPosition());
                    atom.getPosition().PE(originNine);
                    Vector shift = box.getBoundary().centralImage(atom.getPosition());
                    atom.getPosition().PE(shift);
                });
                IMolecule moleculeMOPTen = box.getMoleculeList().get(9);
                Vector originTen = new Vector3D(0,0, -50);
                moleculeMOPTen.getChildList().forEach(atom -> {
                    oldPositions.get(atom.getIndex()).E(atom.getPosition());
                    atom.getPosition().PE(originTen);
                    Vector shift = box.getBoundary().centralImage(atom.getPosition());
                    atom.getPosition().PE(shift);
                });
                IMolecule moleculeMOPEleven = box.getMoleculeList().get(10);
                Vector originEleven = new Vector3D(0,0, -60);
                moleculeMOPEleven.getChildList().forEach(atom -> {
                    oldPositions.get(atom.getIndex()).E(atom.getPosition());
                    atom.getPosition().PE(originEleven);
                    Vector shift = box.getBoundary().centralImage(atom.getPosition());
                    atom.getPosition().PE(shift);
                });
                IMolecule moleculeMOPTwelve = box.getMoleculeList().get(11);
                Vector originTwelve = new Vector3D(0,0, 60);
                moleculeMOPTwelve.getChildList().forEach(atom -> {
                    oldPositions.get(atom.getIndex()).E(atom.getPosition());
                    atom.getPosition().PE(originTwelve);
                    Vector shift = box.getBoundary().centralImage(atom.getPosition());
                    atom.getPosition().PE(shift);
                });*/
               /* //box.addNewMolecule(speciesGrapheneOne);
                box.addNewMolecule(speciesGrapheneTwo);
              //  box.addNewMolecule(speciesGrapheneThree);
                box.addNewMolecule(speciesGrapheneFour);
                box.addNewMolecule(speciesGrapheneFive);
               // box.addNewMolecule(speciesGrapheneSix);
                box.addNewMolecule(speciesGrapheneSeven);
                box.addNewMolecule(speciesGrapheneEight);
                box.addNewMolecule(speciesGrapheneNine);
              //  box.addNewMolecule(speciesGrapheneTen);
              //  box.addNewMolecule(speciesGrapheneEleven);
                box.addNewMolecule(speciesGrapheneTwelve);
                box.addNewMolecule(speciesGrapheneThirteen);*/
            }else {
                box.addNewMolecule(speciesGrapheneOne);
                box.addNewMolecule(speciesGrapheneFive);
            }
        }
        box.addNewMolecule(speciesMOP);
        //Grapahene and SecondGas
      /* if(ifGraphenePresent==true && ifSecondGasPresent==true){
            if(ifMultipleGraphenePresent){
                sm = new SpeciesManager.Builder().addSpecies(speciesMOP).addSpecies(speciesGrapheneOne).addSpecies(speciesGrapheneFive).addSpecies(speciesGrapheneTwo).addSpecies(speciesGrapheneThree).addSpecies(speciesGrapheneFour). addSpecies(speciesGrapheneSix). addSpecies(speciesGrapheneSeven). addSpecies(speciesGrapheneEight).addSpecies(speciesGas).addSpecies(speciesSecondGas).build();
            }else {
                sm = new SpeciesManager.Builder().addSpecies(speciesMOP).addSpecies(speciesGrapheneOne). addSpecies(speciesGrapheneFive).addSpecies(speciesGas).addSpecies(speciesSecondGas).build();
            }
            // sm = new SpeciesManager.Builder().addSpecies(speciesGrapheneOne).addSpecies(speciesGrapheneTwo).addSpecies(speciesGas).build();
        } else if (ifGraphenePresent==true && ifSecondGasPresent==false) {
            if(ifMultipleGraphenePresent){
                sm = new SpeciesManager.Builder().addSpecies(speciesMOP).addSpecies(speciesGrapheneOne).addSpecies(speciesGrapheneTwo).addSpecies(speciesGrapheneThree).addSpecies(speciesGrapheneFour). addSpecies(speciesGrapheneFive). addSpecies(speciesGrapheneSix). addSpecies(speciesGrapheneSeven). addSpecies(speciesGrapheneEight).addSpecies(speciesGas).build();
            }else {
                sm = new SpeciesManager.Builder().addSpecies(speciesMOP).addSpecies(speciesGrapheneOne). addSpecies(speciesGrapheneFive).addSpecies(speciesGas).build();
            }
        } else if (ifGraphenePresent==false && ifSecondGasPresent==true) {
            sm = new SpeciesManager.Builder().addSpecies(speciesMOP).addSpecies(speciesGas).addSpecies(speciesSecondGas).build();
        } else {
            sm = new SpeciesManager.Builder().addSpecies(speciesMOP).addSpecies(speciesGas).build();
        }*/
       // sm = new SpeciesManager.Builder().addSpecies(speciesMOP).addSpecies(speciesGrapheneOne). addSpecies(speciesGrapheneFive).addSpecies(speciesGas).addSpecies(speciesSecondGas).build();
      // sm = new SpeciesManager.Builder().addSpecies(speciesMOP).addSpecies(speciesGas).build();
        //sm = new SpeciesManager.Builder().addSpecies(speciesGrapheneOne).addSpecies(speciesGrapheneFive).addSpecies(speciesGrapheneTwo).addSpecies(speciesGrapheneThree).addSpecies(speciesGrapheneFour). addSpecies(speciesGrapheneSix). addSpecies(speciesGrapheneSeven). addSpecies(speciesGrapheneEight).addSpecies(speciesGrapheneThirteen).addSpecies(speciesGas).addSpecies(speciesSecondGas).build();
       // sm = new SpeciesManager.Builder().addSpecies(speciesGrapheneOne).addSpecies(speciesGrapheneFive).addSpecies(speciesGrapheneTwo).addSpecies(speciesGrapheneThree).addSpecies(speciesGrapheneFour). addSpecies(speciesGrapheneSix). addSpecies(speciesGrapheneSeven). addSpecies(speciesGrapheneEight).addSpecies(speciesGrapheneNine).addSpecies(speciesGrapheneTen).addSpecies(speciesGrapheneEleven).addSpecies(speciesGrapheneTwelve).addSpecies(speciesGrapheneThirteen).addSpecies(speciesGas).addSpecies(speciesSecondGas).build();
        //Bonding parameters MOP
        //sm = new SpeciesManager.Builder().addSpecies(speciesMOP).addSpecies(speciesGrapheneOne).addSpecies(speciesGrapheneFive).addSpecies(speciesGrapheneTwo).addSpecies(speciesGrapheneThree).addSpecies(speciesGrapheneFour). addSpecies(speciesGrapheneSix). addSpecies(speciesGrapheneSeven). addSpecies(speciesGrapheneEight).addSpecies(speciesGas).addSpecies(speciesSecondGas).build();
//.addSpecies(speciesGrapheneFour). addSpecies(speciesGrapheneSeven).addSpecies(speciesGrapheneOne).addSpecies(speciesGrapheneThree).addSpecies(speciesGrapheneSix).addSpecies(speciesGrapheneTen).addSpecies(speciesGrapheneEleven)
        sm = new SpeciesManager.Builder().addSpecies(speciesMOP).addSpecies(speciesGas).build();
      /*  ArrayList<ArrayList<Integer>> connectedAtoms1 = pdbReaderMOP.getConnectivity();
        ArrayList<ArrayList<Integer>> connectivityModified1 = pdbReaderMOP.getConnectivityModified();
        Map<Integer,String> atomMap1 = pdbReaderMOP.getAtomMapWithoutRunning();
        Map<Integer, String> atomMapModified1 = pdbReaderMOP.getAtomMapModified();
        ArrayList<Integer> bondList1 = pdbReaderMOP.getBondList(connectedAtoms1, atomMap1);
        Map<String, double[]> atomicPotMap1 = pdbReaderMOP.atomicPotMap();
        ArrayList<Integer> bondsNum1 = pdbReaderMOP.getBonds();
        Map<Integer, String> atomIdentifierMapModified1 = pdbReaderMOP.getModifiedAtomIdentifierMap();
        List<int[]> dupletsSorted1= pdbReaderMOP.getDupletesSorted();
        List<int[]>tripletsSorted1= pdbReaderMOP.getAnglesSorted();
        List<int[]>quadrupletsSorted1= pdbReaderMOP.getTorsionSorted();
        //Map<String[],List<int[]>> bondTypesMap1= pdbReaderMOP.idenBondTypes(dupletsSorted1, atomIdentifierMapModified1);
       // Map<String[],List<int[]>> angleTypesMap1= pdbReaderMOP.idenAngleTypes(tripletsSorted1, atomIdentifierMapModified1);
       // Map<String[],List<int[]>> torsionTypesMap1= pdbReaderMOP.idenTorsionTypes(quadrupletsSorted1, atomIdentifierMapModified1);*/
       // System.out.println(speciesGas.getMass()/speciesMOP.getMass());

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
        Map<String[],List<int[]>> bondTypesMap2 = null;
        Map<String[],List<int[]>> angleTypesMap2= null;
        Map<String[],List<int[]>> torsionTypesMap2= null;

        //Bonding Parameters GasOne
        if(!isGasCOMPASS && !isGasTraPPE){
            connectedAtoms2 =pdbReaderReplica.getConnectivityWithoutRunning();
            connectivityModified2 = pdbReaderReplica.getConnectivityModifiedWithoutRunning();
            atomMap2 = pdbReaderReplica.getAtomMapWithoutRunning();
            bondList2 = pdbReaderReplica.getBondList(connectedAtoms2, atomMap2);
            atomicPotMap2 = pdbReaderReplica.atomicPotMap();
            atomIdentifierMapModified2 = pdbReaderReplica.getatomIdentifierMapModified();
            dupletsSorted2= pdbReaderReplica.getDupletesSorted();
            tripletsSorted2=pdbReaderReplica.getAnglesSorted();
            quadrupletsSorted2=pdbReaderReplica.getTorsionSorted();
            bondsNum2 = pdbReaderReplica.getBonds();
            bondTypesMap2= pdbReaderReplica.idenBondTypes(dupletsSorted2, atomIdentifierMapModified2);
            angleTypesMap2= pdbReaderReplica.idenAngleTypes(tripletsSorted2, atomIdentifierMapModified2);
            torsionTypesMap2= pdbReaderReplica.idenTorsionTypes(quadrupletsSorted2, atomIdentifierMapModified2);
        } else {
            connectedAtoms2 =pdbReaderCOMPASS.getConnectivity();
           // connectivityModified2 = pdbReaderCOMPASS.getConnectivityModifiedWithoutRunning();
            //atomMap2 = pdbReaderCOMPASS.getAtomMapWithoutRunning();
            //bondList2 = pdbReaderCOMPASS.getBondList(connectedAtoms2, atomMap2);
           // atomicPotMap2 = pdbReaderCOMPASS.atomicPotMap();
            atomIdentifierMapModified2 = pdbReaderCOMPASS.getAtomIdentifierMapMod();
            /*dupletsSorted2= pdbReaderCOMPASS.getDupletesSorted();
            tripletsSorted2=pdbReaderCOMPASS.getAnglesSorted();
            quadrupletsSorted2=pdbReaderCOMPASS.getTorsionSorted();
            bondsNum2 = pdbReaderCOMPASS.getBonds();
            bondTypesMap2= pdbReaderCOMPASS.idenBondTypes(dupletsSorted2, atomIdentifierMapModified2);
            angleTypesMap2= pdbReaderCOMPASS.idenAngleTypes(tripletsSorted2, atomIdentifierMapModified2);
            torsionTypesMap2= pdbReaderCOMPASS.idenTorsionTypes(quadrupletsSorted2, atomIdentifierMapModified2);*/
        }


        //SetBonding Potential
        UniversalSimulation.makeAtomPotentials(sm);
        PotentialMasterBonding pmBonding = new PotentialMasterBonding(sm, box);
        PotentialMasterCell potentialMasterCell = new PotentialMasterCell(getSpeciesManager(), box, 5, pmBonding.getBondingInfo());
        //MOP
        // SetPotential.SetBondStrech(speciesMOP, bondTypesMap1, angleTypesMap1, torsionTypesMap1,bondsNum1,bondList1, quadrupletsSorted1, atomIdentifierMapModified1,atomicPotMap1, pmBonding);
        //GasOne
        if(atomIdentifierMapModified2.size() > 1){
            SetPotential.setBondStretch(speciesGas,bondTypesMap2, angleTypesMap2, torsionTypesMap2, bondsNum2, bondList2,quadrupletsSorted2, atomIdentifierMapModified2, atomicPotMap2, pmBonding);
        }

        //GasTwo
        if(ifSecondGasPresent){
            SetPotential.setBondStretch(speciesSecondGas,bondTypesMapSecondGas, angleTypesMapSecondGas, torsionTypesMapSecondGas, bondsNumSecondGas, bondListSecondGas,quadrupletsSortedSecondGas, atomIdentifierMapModifiedSecondGas, atomicPotMapSecondGas, pmBonding);
        }
        List<AtomType> listGas;

        //NonBonded
        potentialMasterCell.doAllTruncationCorrection = false;
        if(ifSecondGasPresent){
            listGas = SetPotential.listTwoSpeciesPairs(speciesGas.getUniqueAtomTypes(), speciesSecondGas.getUniqueAtomTypes());
        } else {
            listGas = speciesGas.getUniqueAtomTypes();
        }

       // listMOPGasPairs = SetPotential.listGrapheneSpecial(speciesMOP.getUniqueAtomTypes(), speciesGas.getUniqueAtomTypes());
       // List<List<AtomType>> listMOPGasPairs = SetPotential.listFinal(listGas);


      //  List<AtomType> listMOPGas = SetPotential.listTwoSpeciesPairs(speciesMOP.getUniqueAtomTypes(), listGas);
        List<List<AtomType>> listMOPGasPairs = SetPotential.listGrapheneSpecial(speciesMOP.getUniqueAtomTypes(), listGas);
       // List<List<AtomType>> listMOPGasPairs = SetPotential.listFinal( listGas);
        //MOP-Gas

       // System.out.println(listMOPGas);
        LJUFF[] p2LJMOPGas = new LJUFF[listMOPGasPairs.size()];
        P2Electrostatic[] p2ElectroMOPGas = new P2Electrostatic[listMOPGasPairs.size()];
        LJCOMPASS[] p2LJMOPGasCOMPASS = new LJCOMPASS[listMOPGasPairs.size()];
        if(isGasTraPPE){
            SetPotential.doLJElectrostatic(listMOPGasPairs, potentialMasterCell,  p2LJMOPGas, p2ElectroMOPGas, listMOPGasPairs.size(), truncatedRadiusLJ, false, true);
        }else {
            SetPotential.doLJElectrostatic(listMOPGasPairs, potentialMasterCell,  p2LJMOPGas, p2ElectroMOPGas, listMOPGasPairs.size(), truncatedRadiusLJ, false);
        }


        //Graphene-Gas
        if(ifGraphenePresent){
            List<List<AtomType>> listGrapheneGasFinal = SetPotential.listGrapheneSpecial(speciesGrapheneOne.getUniqueAtomTypes(), listGas);
          //  List<List<AtomType>> listGrapheneGasFinal = SetPotential.listFinal(listGrapheneGasUnique);
            LJUFF[] p2LJGrapheneGas = new LJUFF[listGrapheneGasFinal.size()];
            SetPotential.doLJ(listGrapheneGasFinal, potentialMasterCell,  p2LJGrapheneGas, listGrapheneGasFinal.size(), truncatedRadiusLJ);
        }

        integrator = new IntegratorMC(potentialMasterCell, random, temperature, box);
        potentialMasterCell.doOneTruncationCorrection = true;
      //  potentialMasterCell.init();
        if(ifMoveRotateMoves){
            mcMoveMolecule = new MCMoveMolecule(random, potentialMasterCell, box);
            mcMoveMolecule.setStepSize(0.2 * sigma);
            ((MCMoveStepTracker) mcMoveMolecule.getTracker()).setTunable(false);
            mcMoveMoleculeRotate = new MCMoveMoleculeRotate(random, potentialMasterCell, box);
            mcMoveMoleculeRotate.setStepSize(0.2 * sigma);
            ((MCMoveStepTracker) mcMoveMoleculeRotate.getTracker()).setTunable(false);

            if(makeAllMove){
                ((MoleculeSourceRandomMolecule) mcMoveMolecule.getMoleculeSource()).setSpecies(speciesMOP);
                ((MoleculeSourceRandomMolecule) mcMoveMoleculeRotate.getMoleculeSource()).setSpecies(speciesMOP);
            }else {
                ((MoleculeSourceRandomMolecule) mcMoveMolecule.getMoleculeSource()).setSpecies(speciesGas);
               // ((MoleculeSourceRandomMolecule) mcMoveMolecule.getMoleculeSource()).setSpecies(speciesMOP);
               // ((MoleculeSourceRandomMolecule) mcMoveMoleculeRotate.getMoleculeSource()).setSpecies(speciesMOP);
                ((MoleculeSourceRandomMolecule) mcMoveMoleculeRotate.getMoleculeSource()).setSpecies(speciesGas);
            }
            ((MoleculeSourceRandomMolecule) mcMoveMolecule.getMoleculeSource()).setSpecies(speciesGas);
            integrator.getMoveManager().addMCMove(mcMoveMolecule);
            integrator.getMoveManager().addMCMove(mcMoveMoleculeRotate);
        }

        mcMoveID = new MCMoveInsertDelete(potentialMasterCell, random, space);
        mcMoveID.setBox(box);
        mcMoveID.setMu(mu1);
        mcMoveID.setSpecies(speciesGas);
        integrator.getMoveManager().addMCMove(mcMoveID);

        if(ifSecondGasPresent){
            mcMoveIDTwo = new MCMoveInsertDelete(potentialMasterCell, random,space );
            mcMoveIDTwo.setBox(box);
            mcMoveIDTwo.setMu(mu2); //incorrect value
            mcMoveIDTwo.setSpecies(speciesSecondGas);
            integrator.getMoveManager().addMCMove(mcMoveIDTwo);
        }


        final Vector bs = box.getBoundary().getBoxSize();
        System.out.println(bs);
        System.out.println(confNameGasOne);
        integrator.getMoveManager().setEquilibrating(false);

        BoxInflate inflater = new BoxInflate(box, space);
        inflater.actionPerformed();
        potential = new P2LennardJones(sigma, 1.0);
       //
        if (truncatedRadius > 0.5 * box.getBoundary().getBoxSize().getX(1)) {
            throw new RuntimeException("Truncation radius too large.  Max allowed is" + 0.5 * box.getBoundary().getBoxSize().getX(0));
        }
    }

    public static void main(String[] args) {
        GCMCMOPParams params = new GCMCMOPParams();
        ParseArgs.doParseArgs(params, args);
        int numAtomOne = params.numAtomOne;
        //int numAtomTwo = params.numAtomTwo;
        String confNameOne = params.confNameOne;
        String confNameGasOne = params.confNameGasOne;
        String confNameGasTwo = params.confNameGasTwo;
        String confNameGraphene = params.confNameGraphene;
        Vector boxSize = params.boxSize;
        double truncatedRadius = params.truncatedRadius;
        double truncatedRadiusLJ = params.truncatedRadiusLJ;
        double temperature = params.temperature;
        Vector centreMOP = params.centreMOP;
        Vector grapheneOne = params.grapheneOne;
        Vector grapheneTwo = params.grapheneTwo;
        Vector grapheneThree = params.grapheneThree;
        Vector grapheneFour = params.grapheneFour;
        Vector grapheneFive = params.grapheneFive;
        Vector grapheneSix = params.grapheneSix;
        Vector grapheneSeven = params.grapheneSeven;
        Vector grapheneEight = params.grapheneEight;
        Vector grapheneNine = params.grapheneNine;
        Vector grapheneTen = params.grapheneTen;
        Vector grapheneEleven = params.grapheneEleven;
        Vector grapheneTwelve = params.grapheneTwelve;
        Vector grapheneThirteen = params.grapheneThirteen;
        Vector centreMOPTwo = params.centreMOPTwo;
        Vector centreMOPThree = params.centreMOPThree;
        Vector centreMOPFour = params.centreMOPFour;
        boolean ifGraphenePresent = params.ifGraphenePresent;
        boolean ifSecondGasPresent = params.ifSecondGasPresent;
        boolean ifMultipleGraphenePresent = params.ifMultipleGraphenePresent;
        boolean ifMoveRotateMoves = params.ifMoveRotateMoves;
        boolean makeAllMove = params.makeAllMove;
        boolean doGraphics = params.doGraphics;
        boolean doElectrostatics = params.doElectrostatics;
       // double mu1 = -params.mu1*BOLTZMANN_K*Kelvin.UNIT.toSim(temperature);
        boolean isGasCOMPASS = params.isGasCOMPASS;
        boolean isGasTraPPE = params.isGasTraPPE;
        double mu1 = params.mu1;
        double multiplier = params.multiplier;
        //mu = -KbTexp(-BU)
        double mu2= params.mu2;
        double sigma = params.sigma;
        System.out.println(mu1 + " mu1");

        GCMCMOP sim = new GCMCMOP(confNameOne, confNameGasOne, confNameGasTwo, confNameGraphene, centreMOP,centreMOPTwo, centreMOPThree, centreMOPFour, grapheneOne, grapheneTwo, grapheneThree, grapheneFour, grapheneFive, grapheneSix, grapheneSeven, grapheneEight,grapheneNine, grapheneTen, grapheneEleven, grapheneTwelve, grapheneThirteen,numAtomOne,  temperature, truncatedRadius, truncatedRadiusLJ,sigma, mu1,mu2, ifGraphenePresent, ifSecondGasPresent,ifMultipleGraphenePresent,ifMoveRotateMoves, boxSize, makeAllMove, doElectrostatics, multiplier, isGasCOMPASS, isGasTraPPE);
        if(doGraphics){
            final SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, "GCMC");
            DiameterHashByType dhbt = (DiameterHashByType) simGraphic.getDisplayBox(sim.box).getDiameterHash();
           /* for(int i=0; i<sim.speciesGrapheneOne.getUniqueAtomTypes().size(); i++){
                if(i==1 || i==2 || i ==3 || i==4){
                    dhbt.setDiameter(sim.speciesGrapheneOne.getAtomType(i), 1.5);
                }else{
                    dhbt.setDiameter(sim.speciesGrapheneOne.getAtomType(i), 1);

                }
            }*/

           // dhbt.setDiameter(sim.speciesMOP.getAtomType(0), 1);
            //dhbt.setDiameter(sim.speciesMOP.getAtomType(1), 1);
            //dhbt.setDiameter(sim.speciesMOP.getAtomType(2), 1);
           // dhbt.setDiameter(sim.speciesMOP.getAtomType(3), 1);
            //dhbt.setDiameter(sim.speciesMOP.getAtomType(4), 1);
            //dhbt.setDiameter(sim.speciesMOP.getAtomType(4), 0.2);
            //dhbt.setDiameter(sim.speciesMOP.getAtomType(5), 0.2);
            if(ifGraphenePresent) {
                /*for (int i=0; i<sim.speciesGrapheneOne.getUniqueAtomTypes().size(); i++){
                    dhbt.setDiameter(sim.speciesGrapheneOne.getAtomType(i), 1.5);
                }*/
               /* dhbt.setDiameter(sim.speciesGas.getAtomType(0), 1.5);
                dhbt.setDiameter(sim.speciesGas.getAtomType(1), 1);
                dhbt.setDiameter(sim.speciesGrapheneOne.getAtomType(0), 1.2);
              dhbt.setDiameter(sim.speciesGrapheneOne.getAtomType(1), 1.2);
                dhbt.setDiameter(sim.speciesGrapheneOne.getAtomType(2), 1.2);
                dhbt.setDiameter(sim.speciesGrapheneOne.getAtomType(3), 0.8);
                dhbt.setDiameter(sim.speciesGrapheneOne.getAtomType(4), 1.2);
                dhbt.setDiameter(sim.speciesGrapheneOne.getAtomType(5), 1.2);
                dhbt.setDiameter(sim.speciesGrapheneOne.getAtomType(6), 1.2);
                dhbt.setDiameter(sim.speciesGrapheneOne.getAtomType(7), 1.2);
                dhbt.setDiameter(sim.speciesGrapheneOne.getAtomType(8), 1.2);
                ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGrapheneOne.getTypeByName("CX"), Color.darkGray);
                ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGrapheneOne.getTypeByName("CY"), Color.darkGray);
                ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGrapheneOne.getTypeByName("CZ"), Color.darkGray);
                ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGrapheneOne.getTypeByName("C4"), ColorExtra.darkGray);
                ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGrapheneOne.getTypeByName("OL"), ColorExtra.oliveDrab);
                ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGrapheneOne.getTypeByName("OJ"), ColorExtra.oliveDrab);
                ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGrapheneOne.getTypeByName("OK"), ColorExtra.oliveDrab);
                ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGrapheneOne.getTypeByName("OE"), ColorExtra.oliveDrab);
                ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGrapheneOne.getTypeByName("HK"), ColorExtra.lightcyan1);
                ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGrapheneOne.getTypeByName("C4"), Color.darkGray);
               // ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGrapheneOne.getTypeByName("HK"), ColorExtra.Violet);
                //((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGrapheneOne.getTypeByName("OJ"), ColorExtra.gold);
               // ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGrapheneOne.getTypeByName("OK"), ColorExtra.cornflowerblue);
          /*      for (int i=0; i<sim.speciesMOP.getUniqueAtomTypes().size();i++){
                    ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesMOP.getAtomType(i), ColorExtra.white);
                }


            }
            System.out.println(sim.speciesMOP.getMass());*/
            }
            if (ifGraphenePresent){
                ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGrapheneOne.getAtomType(0), Color.gray);
              ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGrapheneOne.getAtomType(1), Color.gray);
                ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGrapheneOne.getAtomType(2), ColorExtra.lightcyan);
                ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGrapheneOne.getAtomType(3), Color.red);
                ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGrapheneOne.getAtomType(4), Color.gray);
                ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGrapheneOne.getAtomType(5), ColorExtra.lightcyan);
                ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGrapheneOne.getAtomType(6), ColorExtra.lightcyan);
                ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGrapheneOne.getAtomType(7), ColorExtra.lightcyan);
            }
           // dhbt.setDiameter(sim.speciesGas.getAtomType(1),1);


            //MOP
        /*   dhbt.setDiameter(sim.speciesMOP.getAtomType(0), 3);
            dhbt.setDiameter(sim.speciesMOP.getAtomType(1), 2);
            dhbt.setDiameter(sim.speciesMOP.getAtomType(2), 1);
            dhbt.setDiameter(sim.speciesMOP.getAtomType(3), 2);
            dhbt.setDiameter(sim.speciesMOP.getAtomType(4), 2);*/
           // dhbt.setDiameter(sim.speciesMOP.getAtomType(5), 2);
           // dhbt.setDiameter(sim.speciesMOP.getAtomType(5), 2);
           // dhbt.setDiameter(sim.speciesGas.getAtomType(0), 1.5);
           // dhbt.setDiameter(sim.speciesGas.getAtomType(1), 1);
          //  dhbt.setDiameter(sim.speciesSecondGas.getAtomType(0), 1.5);
           // dhbt.setDiameter(sim.speciesSecondGas.getAtomType(1), 1);
          //  ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesMOP.getAtomType(0), ColorExtra.copper);
           // ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesMOP.getAtomType(1), ColorExtra.red);
            //((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesMOP.getAtomType(2), ColorExtra.blue);
           // ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesMOP.getAtomType(3), ColorExtra.gray);
            //((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesMOP.getAtomType(4), ColorExtra.firebrick);
          // ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesMOP.getAtomType(5), ColorExtra.gray);

            //((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGas.getAtomType(0), ColorExtra.lightseagreen);
           // ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGas.getAtomType(1), ColorExtra.salmon);

        /*   ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesMOP.getTypeByName("Cu"), Color.white);
            ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesMOP.getTypeByName("H"), Color.RED);
          ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesMOP.getTypeByName("C_3"), Color.cyan);
            ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesMOP.getTypeByName("O_3"), Color.blue);
            ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesMOP.getTypeByName("O_1"), Color.blue);*/
          // ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesMOP.getAtomType(0), ColorExtra.lightcyan);
           /*((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGas.getAtomType(1), ColorExtra.salmon);
            ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesSecondGas.getAtomType(0), ColorExtra.cyan);
            ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesSecondGas.getAtomType(1), ColorExtra.firebrick);*/

            simGraphic.makeAndDisplayFrame("GCMC");
            ActivityIntegrate ai2 = new ActivityIntegrate(sim.integrator);
            sim.getController().addActivity(ai2, Long.MAX_VALUE, 1.0);
            return;
        }
        //System.out.println(sim.box.getMoleculeList(sim.speciesGrapheneOne).size()* sim.speciesGrapheneOne.getMass()+ " " + sim.box.getMoleculeList(sim.speciesGrapheneOne).size()+ " " + sim.speciesGrapheneOne.getMass());

        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, params.numSteps / 10));

        int pInterval = 2 ;
        int bs = params.numSteps / (pInterval * 5);
        if (bs == 0) bs = 1;
        MeterPressure pMeter = new MeterPressure(sim.box, sim.integrator.getPotentialCompute());
        pMeter.setTemperature(sim.integrator.getTemperature());
        AccumulatorAverage pAccumulator = new AccumulatorAverageFixed(bs);
        DataPumpListener pPump = new DataPumpListener(pMeter, pAccumulator, pInterval);
        sim.integrator.getEventManager().addListener(pPump);

        bs = params.numSteps / 50;
        if (bs == 0) bs = 1;
        MeterPotentialEnergyFromIntegrator energyMeter = new MeterPotentialEnergyFromIntegrator(sim.integrator);
        AccumulatorAverage energyAccumulator = new AccumulatorAverageFixed(bs);
        DataPumpListener uPump = new DataPumpListener(energyMeter, energyAccumulator);
        sim.integrator.getEventManager().addListener(uPump);

        nMoleculesGasOne = new MeterNMolecules();
        nMoleculesGasOne.setBox(sim.box);
        nMoleculesGasOne.setSpecies(sim.speciesGas);

        if(ifSecondGasPresent){
            nMoleculesGasTwo = new MeterNMolecules();
            nMoleculesGasTwo.setBox(sim.box);
            nMoleculesGasTwo.setSpecies(sim.speciesSecondGas);
        }

        MeterDensity densityMeter = new MeterDensity(sim.box);
        AccumulatorAverage densityAccumulator = new AccumulatorAverageFixed(bs);
        densityMeter.setSpecies(sim.speciesGas);
        DataPumpListener pumpDensity = new DataPumpListener(densityMeter, densityAccumulator);
        sim.integrator.getEventManager().addListener(pumpDensity);

        if(ifSecondGasPresent){
            MeterDensity densityMeterGasTwo = new MeterDensity(sim.box);
            AccumulatorAverage densityAccumulatorGastwo = new AccumulatorAverageFixed(bs);
            densityMeterGasTwo.setSpecies(sim.speciesSecondGas);
            DataPumpListener pumpDensityGasTwo = new DataPumpListener(densityMeterGasTwo, densityAccumulatorGastwo);
            sim.integrator.getEventManager().addListener(pumpDensityGasTwo);
        }


        long t1 = System.currentTimeMillis();
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, params.numSteps));
        long t2 = System.currentTimeMillis();
        System.out.println("mu1 " + mu1);
        System.out.println("runtime: " + (t2 - t1) * 0.001);
        Unit pUnit = Bar.UNIT;
        Unit MegaPascal = new PrefixedUnit(Prefix.MEGA, Pascal.UNIT);
        double avgP = pAccumulator.getData(pAccumulator.AVERAGE).getValue(0);
        double errP = pAccumulator.getData(pAccumulator.ERROR).getValue(0);
        double corP = pAccumulator.getData(pAccumulator.BLOCK_CORRELATION).getValue(0);
        System.out.println("P " + MegaPascal.fromSim(avgP) + " " + MegaPascal.fromSim(errP) + " " + corP);

      //  double muR = Constants.BOLTZMANN_K * temperature * Math.log(avgP/temperature);
        //System.out.println("muR: " + muR);

        double avgRho = densityAccumulator.getData(densityAccumulator.AVERAGE).getValue(0);
        double errRho = densityAccumulator.getData(densityAccumulator.ERROR).getValue(0);
        double corRho = densityAccumulator.getData(densityAccumulator.BLOCK_CORRELATION).getValue(0);
        Unit dm = new PrefixedUnit(Prefix.DECI, Meter.UNIT);
        Unit dm3 = new CompoundUnit(new Unit[] {dm}, new double[] {3});
        Unit m3 = new CompoundUnit(new Unit[] {Meter.UNIT}, new double[] {3});
        Unit kgm3 = new UnitRatio(new PrefixedUnit(Prefix.KILO, Gram.UNIT), m3);
        Unit moldm3 = new UnitRatio(Mole.UNIT, Liter.UNIT);
        double numAtomsAvg = avgRho * sim.box.getBoundary().volume();
        double errnum = errRho *  sim.box.getBoundary().volume();
        double cornum = corRho;
        System.out.println( "Num atoms : " +numAtomsAvg + " " + errnum + " " +cornum );
        //System.out.println("Rho (mol/dm3) "+kgm3.fromSim(avgRho)/sim.speciesGas.getMass() + " " + sim.speciesGas.getMass());
        AtomType C1 = new AtomType(Carbon.INSTANCE, "C1");
        System.out.println("Amount sorbed " + numAtomsAvg*100*sim.speciesGas.getMass() /(sim.speciesMOP.getMass()) + " err " +errnum*100*sim.speciesGas.getMass() /(sim.speciesMOP.getMass()));
       // System.out.println("Amount sorbed " + 1000*numAtomsAvg*sim.speciesGas.getMass()/sim.speciesMOP.getMass()  + " err " + 1000*errnum*sim.speciesGas.getMass()/sim.speciesMOP.getMass());
        // System.out.println("Mass fraction : "+sim.speciesMOP.getMass()*100/(sim.speciesGrapheneSix.getMass()*12));
      //  double avgPE = energyAccumulator.getData(energyAccumulator.AVERAGE).getValue(0) / numAtomsAvg;
      //  double errPE = energyAccumulator.getData(energyAccumulator.ERROR).getValue(0) / numAtomsAvg;
      //  double corPE = energyAccumulator.getData(energyAccumulator.BLOCK_CORRELATION).geValue(0);
        //System.out.println("PE " + avgPE + " " + errPE + " " + corPE);


        //double numAtomGasOne = nMoleculesGasOne.getDataAsScalar();
        //System.out.println("Atoms gas One: " +numAtomGasOne);

        if(ifSecondGasPresent){
            double numAtomGasTwo = nMoleculesGasTwo.getDataAsScalar();
            System.out.println("Atoms gas Two: " +numAtomGasTwo);
        }



        System.out.println("rho " + kgm3.fromSim(avgRho) + " " + errRho + " " + corRho);
        System.out.println("uptake  " + 22.414*kgm3.fromSim(avgRho)/sim.speciesGas.getMass());
        //System.out.println("Rho (mol/dm3) " + moldm3.fromSim(avgRho) + " " +errRho + " " +corRho);
        // expected values based on 10^8 steps
        // stdev based on 10^8 steps for N=500 uncertainty, scaled up to 10^6 steps
        //   other sims have significant correlation
        // 4 sigma should fail 1 in 16,000 runs

        double expectedP = 0.0815; // finite size effect smaller than uncertainty
        double stdevP = 0.03;
        if (Double.isNaN(avgP) || Math.abs(avgP - expectedP) / stdevP > 4) {
            System.exit(1);
        }

        //double expectedPE = -4.492; // finite size effect comparable to uncertainty
        //double stdevPE = 0.09;
       /* if (Double.isNaN(avgPE) || Math.abs(avgPE - expectedPE) / stdevPE > 4) {
            System.exit(2);
        }*/

        double expectedRho = 0.65; // finite size effect smaller than uncertainty
        double stdevRho = 0.006;
        if (Double.isNaN(avgRho) || Math.abs(avgRho - expectedRho) / stdevRho > 4) {
            System.exit(3);
        }
    }
   public Integrator getIntegrator(){
       return integrator;
   }
    public static class GCMCMOPParams extends ParameterBase {
      //  public Vector centreMOP = new Vector3D(0.0,0.0, 0.0);
        public Vector boxSize = new Vector3D(30,30,30);
        public Vector grapheneThirteen = new Vector3D(15.0,15.0, 60.0);
        public Vector grapheneOne = new Vector3D(15.0,15.0, 0.0);
        public Vector grapheneTwo = new Vector3D(0,0, 12.0);
        public Vector grapheneThree = new Vector3D(15.0,15.0, 20.0);
        public Vector grapheneFour = new Vector3D(15.0,15.0, 30.0);
        public Vector grapheneNine = new Vector3D(15.0,15.0, 40.0);
        public Vector grapheneTen = new Vector3D(15.0,15.0, 50.0);
        public Vector grapheneFive = new Vector3D(15.0,15.0, -10.0);
        public Vector grapheneSix = new Vector3D(0,0, -12.0);
        public Vector grapheneSeven = new Vector3D(15.0,15.0, -30.0);
        public Vector grapheneEight = new Vector3D(15.0,15.0, -40.0);
        public Vector grapheneEleven = new Vector3D(15.0,15.0, -50.0);
        public Vector grapheneTwelve = new Vector3D(15.0,15.0, -60.0);
        public double mu2 = -2100;
        public Vector centreMOPTwo = new Vector3D(0.0,0.0,30);
        public Vector centreMOP = new Vector3D(0.0,0.0,-30);
        public Vector centreMOPThree = new Vector3D(0.0,0.0,-50);
        public Vector centreMOPFour = new Vector3D(0.0,0.0,50);
        public String confNameGasTwo ="F://Avagadro//molecule//propane";
        public boolean ifGraphenePresent = false;
        public boolean ifMultipleGraphenePresent = true;
        public boolean ifSecondGasPresent = false;
        public boolean ifMoveRotateMoves = true;
        public boolean doElectrostatics = true;
        public int numAtomOne = 1;
        public int truncatedRadiusLJ = 14;
        public double truncatedRadius = 14;

        public double sigma = 5.2;
        public int numSteps = 10000;
        public boolean makeAllMove = false;
        public double mu1 = -2300;
        public double partial = -82.2479;

        public double temperature = 330;
        public String confNameGasOne = "F://Avagadro//molecule//propene" ;
        public String confNameGraphene = "F://Avagadro//holeGO60";
        public String confNameOne = "F://Avagadro//mop//tetra_cu_3C";
        public boolean doGraphics =true;
        public boolean isGasCOMPASS = false;
        public boolean isGasTraPPE = true;
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
