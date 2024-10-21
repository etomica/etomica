package etomica.GasMOP;

import etomica.action.BoxInflate;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.DiameterHashByType;
import etomica.box.Box;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.AccumulatorHistogram;
import etomica.data.DataPumpListener;
import etomica.data.histogram.Histogram;
import etomica.data.histogram.HistogramSimple;
import etomica.data.histogram.HistogramVectorSimple;
import etomica.data.meter.MeterDensity;
import etomica.data.meter.MeterNMolecules;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.meter.MeterPressure;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.Integrator;
import etomica.integrator.IntegratorEvent;
import etomica.integrator.IntegratorListener;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.*;
import etomica.math.DoubleRange;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.molecule.MoleculeSourceRandomMolecule;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.potential.*;
import etomica.potential.TraPPE.SpeciesGasTraPPE;
import etomica.potential.UFF.*;
import etomica.simulation.Simulation;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesGeneral;
import etomica.species.SpeciesManager;
import etomica.units.*;
import etomica.util.Constants;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import java.awt.*;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.List;

import static etomica.potential.TraPPE.SpeciesGasTraPPE.ChemForm.*;

public class GCMCMOP extends Simulation {
    public IntegratorMC integrator;
    public MCMoveMolecule mcMoveMolecule;
    public PotentialMasterCell potentialMasterCell;
    public MCMoveMoleculeRotate mcMoveMoleculeRotate;
    public MCMoveInsertDelete mcMoveID, mcMoveIDTwo;
    public ISpecies speciesMOP, speciesGas, speciesGrapheneOne, speciesGrapheneTwo, speciesGrapheneThree, speciesGrapheneFour,speciesGrapheneFive,speciesGrapheneSix,speciesGrapheneSeven,speciesGrapheneEight, speciesSecondGas,speciesGrapheneNine,speciesGrapheneTen,speciesGrapheneEleven,speciesGrapheneTwelve, speciesGrapheneThirteen;
    public ISpecies speciesMOPTwo, speciesMOPThree,speciesMOPFour, speciesMOPFive, speciesCOF;
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
    public static double massMembrane;
    public static HistogramVectorSimple targHist;
    public GCMCMOP(String confNameOne, String confNameGasOne, String confNameGasTwo, String confNameGraphene, Vector centreMOP, Vector centreMOPTwo, Vector centreMOPThree,  Vector centreMOPFour, Vector grapheneOne, Vector grapheneTwo,Vector grapheneThree,Vector grapheneFour,Vector grapheneFive,Vector grapheneSix,Vector grapheneSeven,Vector grapheneEight,Vector grapheneNine, Vector grapheneTen, Vector grapheneEleven, Vector grapheneTwelve, Vector grapheneThirteen, int numMolOne, int numMolTwo, double temperature, double truncatedRadius, double truncatedRadiusLJ, double sigma, double mu1, double mu2, boolean ifGraphenePresent, boolean ifSecondGasPresent, boolean ifMultipleGraphenePresent, boolean ifMoveRotateMoves, Vector boxSize, boolean makeAllMove, boolean doElectrostatics, double multiplier, boolean isGasCOMPASS, boolean isGasTraPPE, boolean ifMOPPresent, boolean ifCOFPesent, String automMOPGeom, boolean ifautoMOP) {
        super(Space3D.getInstance());

        //Make Species
       // setRandom(new RandomMersenneTwister(2));
        PDBReaderMOP pdbReaderMOP = new PDBReaderMOP();
        PDBReaderMOP pdbReaderMOP2 = new PDBReaderMOP();
        PDBReaderMOP pdbReaderCOF = new PDBReaderMOP();
        PDBReaderReplica pdbReaderReplica = new PDBReaderReplica();
        //PDBReaderCOMPASS pdbReaderCOMPASS = new PDBReaderCOMPASS();
        PDBReaderReplicaNew pdbReaderReplicaNew = new PDBReaderReplicaNew();
        GeneralGrapheneReader grapheneReader = new GeneralGrapheneReader();
        SetPotential SetPotential = new SetPotential();


        if (ifGraphenePresent) {
            if (ifMultipleGraphenePresent) {
                speciesGrapheneOne = grapheneReader.getSpecies(confNameGraphene, grapheneOne, false);
                speciesGrapheneFive = grapheneReader.getSpecies(confNameGraphene, grapheneFive, false);
                addSpecies(speciesGrapheneOne);
                addSpecies(speciesGrapheneFive);
            } else {
                speciesGrapheneOne = grapheneReader.getSpecies(confNameGraphene, grapheneOne, false);
                addSpecies(speciesGrapheneOne);
            }
        }
        if (ifMOPPresent) {
            speciesMOP = pdbReaderMOP.getSpeciesMOP(confNameOne, false, new Vector3D(0.0, 0.0, 0.0), false);
            cifReader cifReader = new cifReader();
            addSpecies(speciesMOP);
        }
      /*  if(ifCOFPesent){
            speciesCOF = pdbReaderCOF.getSpecies(c)
        }*/

        SpeciesGasTraPPE speciesGasTraPPE = new SpeciesGasTraPPE();
        if (isGasTraPPE) {
            if (confNameGasOne.equals("F://Avagadro//molecule//ch4")) {
                SpeciesGasTraPPE.ChemForm = CH4;
            } else if (confNameGasOne.equals("F://Avagadro//molecule//ethane")) {
                SpeciesGasTraPPE.ChemForm = C2H6;
            } else if (confNameGasOne.equals("F://Avagadro//molecule//propane")) {
                SpeciesGasTraPPE.ChemForm = C3H8;
            } else if (confNameGasOne.equals("F://Avagadro//molecule//ethene")) {
                SpeciesGasTraPPE.ChemForm = C2H4;
            } else if (confNameGasOne.equals("F://Avagadro//molecule//propene")) {
                SpeciesGasTraPPE.ChemForm = C3H6;
            } else if (confNameGasOne.equals("F://Avagadro//molecule//co2")) {
                SpeciesGasTraPPE.ChemForm = CO2;
            } else if (confNameGasOne.equals("F://Avagadro//molecule//o2")) {
                SpeciesGasTraPPE.ChemForm = O2;
            } else if (confNameGasOne.equals("F://Avagadro//molecule//n2")) {
                SpeciesGasTraPPE.ChemForm = N2;
            }else if (confNameGasOne.equals("F://Avagadro//molecule//nh3")) {
                SpeciesGasTraPPE.ChemForm = NH3;
            }
        }
        //speciesMOP = cifReader.speciesCIF(confNameOne,false);

        //species Gas
        if (isGasCOMPASS) {
            //  speciesGas = pdbReaderCOMPASS.getSpecies(confNameGasOne, false, true, new Vector3D(0,0,0));
        } else if (isGasTraPPE) {
            speciesGas = speciesGasTraPPE.speciesGasTraPPE(Space3D.getInstance(), SpeciesGasTraPPE.ChemForm, false);
        } else {
            speciesGas = pdbReaderReplica.getSpecies(confNameGasOne, true, new Vector3D(0, 0, 0), false);
        }
        addSpecies(speciesGas);
        //   System.out.println(speciesMOP.getMass());
        //    System.out.println(speciesGas.getMass());

        //SecondGas and Bonding Parameters of Second Gas
      if(ifSecondGasPresent){
            speciesSecondGas = pdbReaderReplicaNew.getSpecies(confNameGasTwo);
            addSpecies(speciesSecondGas);
            pairsAtomsSecondGas = SetPotential.getSpeciesPairs(speciesSecondGas);
            connectedAtomsSecondGas = pdbReaderReplicaNew.getConnectivityWithoutRunning();
            connectivityModifiedSecondGas = pdbReaderReplicaNew.getConnectivityModifiedWithoutRunning();
            atomMapSecondGas = pdbReaderReplicaNew.getAtomMapWithoutRunning();
            bondListSecondGas = pdbReaderReplicaNew.getBondList(connectedAtomsSecondGas, atomMapSecondGas);
            atomicPotMapSecondGas = pdbReaderReplicaNew.atomicPotMap();
            atomIdentifierMapModifiedSecondGas = pdbReaderReplicaNew.getatomIdentifierMapModified();
            dupletsSortedSecondGas = pdbReaderReplicaNew.getDupletesSorted();
            tripletsSortedSecondGas = pdbReaderReplicaNew.getAnglesSorted();
            quadrupletsSortedSecondGas = pdbReaderReplicaNew.getTorsionSorted();
            bondsNumSecondGas = pdbReaderReplicaNew.getBonds();
            bondTypesMapSecondGas= pdbReaderReplicaNew.idenBondTypes(dupletsSortedSecondGas, atomIdentifierMapModifiedSecondGas);
            angleTypesMapSecondGas= pdbReaderReplicaNew.idenAngleTypes(tripletsSortedSecondGas, atomIdentifierMapModifiedSecondGas);
            torsionTypesMapSecondGas= pdbReaderReplicaNew.idenTorsionTypes(quadrupletsSortedSecondGas, atomIdentifierMapModifiedSecondGas);
        }

      box = this.makeBox();
      double truncCal =0;
      if(ifMOPPresent){
          Map<Integer, Vector> listPositions = pdbReaderMOP.getPositions();
          distCalc distCalc = new distCalc();
          Vector boxSizeCalc = distCalc.boxSizeVector(listPositions);
        //  System.out.println(boxSizeCalc + " " + speciesMOP.getMass());
           truncCal = distCalc.truncRad(boxSizeCalc);
          //System.out.println(boxSizeCalc + " " + truncCal);
          box.getBoundary().setBoxSize(new Vector3D(50,50,50));
      } else {
          box.getBoundary().setBoxSize(boxSize);
      }
      Map<Double, List<Integer[]> >distMap = new HashMap<>();
      //tetra= 6, cube = 12 octahedron = 12, dodecahedron = 30, icosahedron = 30;
      if(ifautoMOP){
          List<Vector> oldPositions = new ArrayList<>();
          List<Double> slopeXY = new ArrayList<>();
          List<Double> angleZ = new ArrayList<>();
          getPolyhedra(distMap,slopeXY, angleZ, automMOPGeom );
          if(automMOPGeom.equals("tetra")){
              box.setNMolecules(speciesMOP, 2);
              //second Strand
              IMolecule moleculeMOPZero = box.getMoleculeList().get(0);
              IMolecule moleculeMOPOne = box.getMoleculeList().get(1);
              IMolecule moleculeMOPTwo = box.getMoleculeList().get(2);
              IMolecule moleculeMOPThree = box.getMoleculeList().get(3);
              IMolecule moleculeMOPFour = box.getMoleculeList().get(4);
              IMolecule moleculeMOPFive = box.getMoleculeList().get(5);
              while (oldPositions.size() < moleculeMOPZero.getChildList().size()) {
                  oldPositions.add(space.makeVector());
              }
              int num = leftSideMost(moleculeMOPOne);
              Vector diffOne = new Vector3D(30,30,30);
              System.out.println(moleculeMOPZero.getChildList().get(num).getPosition() + " "+ moleculeMOPOne.getChildList().get(moleculeMOPOne.getChildList().size()-1).getPosition());
              diffOne.Ev1Mv2(moleculeMOPZero.getChildList().get(num).getPosition(), moleculeMOPOne.getChildList().get(moleculeMOPOne.getChildList().size()-1).getPosition());
              System.out.println(diffOne);
              moleculeMOPOne.getChildList().forEach(atom -> {
                  oldPositions.get(atom.getIndex()).E(atom.getPosition());
                  atom.getPosition().PE(diffOne);
               //   System.out.println(diffOne);
                  Vector shift = box.getBoundary().centralImage(atom.getPosition());
                  atom.getPosition().PE(shift);
              });
              for(int u=0; u<moleculeMOPOne.getChildList().size(); u++){
                  System.out.println(moleculeMOPOne.getChildList().get(u).getPosition());
              }

          //    System.exit(1);
          } else if (automMOPGeom.equals("rhombic")) {
              box.setNMolecules(speciesMOP, 12);
          }else if (automMOPGeom.equals("cubic")) {
              box.setNMolecules(speciesMOP, 12);
          }else if (automMOPGeom.equals("octa")) {
              box.setNMolecules(speciesMOP, 30);
          }else if (automMOPGeom.equals("icosa")) {
              box.setNMolecules(speciesMOP, 30);
          }else if (automMOPGeom.equals("dodeca")) {

          }

      }

     /* List<Vector> oldPositions = new ArrayList<>();
      List<Vector> oldPositionsTwo = new ArrayList<>();
      box.setNMolecules(speciesMOP, 4);
        IMolecule moleculeMOPOne = box.getMoleculeList().get(0);
      while (oldPositions.size() < moleculeMOPOne.getChildList().size()) {
          oldPositions.add(space.makeVector());
      }
        Vector originTwo = new Vector3D(0,0,0);
        IMolecule moleculeMOPTwo = box.getMoleculeList().get(1);
      moleculeMOPTwo.getChildList().forEach(atom -> {
          oldPositions.get(atom.getIndex()).E(atom.getPosition());
          atom.getPosition().PE(originTwo);
          Vector shift = box.getBoundary().centralImage(atom.getPosition());
          atom.getPosition().PE(shift);
      });
      Vector originOne = new Vector3D(0,0, 0);
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
        Vector originFour = new Vector3D(0,0,0 );
        moleculeMOPFour.getChildList().forEach(atom -> {
            oldPositions.get(atom.getIndex()).E(atom.getPosition());
            atom.getPosition().PE(originFour);
            Vector shift = box.getBoundary().centralImage(atom.getPosition());
            atom.getPosition().PE(shift);
        });*/
     /*   IMolecule moleculeMOPFive = box.getMoleculeList().get(4);
        Vector originFive = new Vector3D(0,0,-50 );
        moleculeMOPFive.getChildList().forEach(atom -> {
            oldPositions.get(atom.getIndex()).E(atom.getPosition());
            atom.getPosition().PE(originFive);
            Vector shift = box.getBoundary().centralImage(atom.getPosition());
            atom.getPosition().PE(shift);
        });*/
        if(ifGraphenePresent ){
            if(!ifMultipleGraphenePresent ){
                List<Vector> oldPositions = new ArrayList<>();
                massMembrane = speciesGrapheneOne.getMass()*numMolOne;
                box.setNMolecules(speciesGrapheneOne, 2);
                IMolecule moleculeMOPZero = box.getMoleculeList().get(0);
                while (oldPositions.size() < moleculeMOPZero.getChildList().size()) {
                    oldPositions.add(space.makeVector());
                }
                moleculeMOPZero.getChildList().forEach(atom -> {
                    oldPositions.get(atom.getIndex()).E(atom.getPosition());
                    atom.getPosition().PE(grapheneTwo);
                    Vector shift = box.getBoundary().centralImage(atom.getPosition());
                    atom.getPosition().PE(shift);
                });
                IMolecule moleculeMOPOne = box.getMoleculeList().get(1);
                moleculeMOPOne.getChildList().forEach(atom -> {
                    oldPositions.get(atom.getIndex()).E(atom.getPosition());
                    atom.getPosition().PE(grapheneSix);
                    Vector shift = box.getBoundary().centralImage(atom.getPosition());
                    atom.getPosition().PE(shift);
                });
             /*   IMolecule moleculeMOPThree = box.getMoleculeList().get(2);
                moleculeMOPThree.getChildList().forEach(atom -> {
                    oldPositions.get(atom.getIndex()).E(atom.getPosition());
                    atom.getPosition().PE(grapheneThree);
                    Vector shift = box.getBoundary().centralImage(atom.getPosition());
                    atom.getPosition().PE(shift);
                });
                IMolecule moleculeMOPFour = box.getMoleculeList().get(3);
                moleculeMOPFour.getChildList().forEach(atom -> {
                    oldPositions.get(atom.getIndex()).E(atom.getPosition());
                    atom.getPosition().PE(grapheneSeven);
                    Vector shift = box.getBoundary().centralImage(atom.getPosition());
                    atom.getPosition().PE(shift);
                });

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
                });*/
             /*   IMolecule moleculeMOPThree = box.getMoleculeList().get(2);
                moleculeMOPThree.getChildList().forEach(atom -> {
                    oldPositions.get(atom.getIndex()).E(atom.getPosition());
                    atom.getPosition().PE(grapheneThree);
                    Vector shift = box.getBoundary().centralImage(atom.getPosition());
                    atom.getPosition().PE(shift);
                });
                IMolecule moleculeMOPEight = box.getMoleculeList().get(3);
                moleculeMOPEight.getChildList().forEach(atom -> {
                    oldPositions.get(atom.getIndex()).E(atom.getPosition());
                    atom.getPosition().PE(grapheneSeven);
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
            } else {
                box.addNewMolecule(speciesGrapheneOne);
                box.addNewMolecule(speciesGrapheneFive);
            }
        }
        //System.out.println(speciesGas.getMass());
        if(ifMOPPresent ){
            massMembrane = speciesMOP.getMass()*numMolTwo;
            box.setNMolecules(speciesMOP, numMolTwo);
           /* IMolecule moleculeMOPFour = box.getMoleculeList().get(4);
            while (oldPositionsTwo.size() < moleculeMOPFour.getChildList().size()) {
                oldPositionsTwo.add(space.makeVector());
            }

            IMolecule moleculeMOPFive = box.getMoleculeList().get(4);
            Vector originTwelve = new Vector3D(0,-10, 0);
            moleculeMOPFive.getChildList().forEach(atom -> {
                oldPositionsTwo.get(atom.getIndex()).E(atom.getPosition());
                atom.getPosition().PE(originTwelve);
                Vector shift = box.getBoundary().centralImage(atom.getPosition());
                atom.getPosition().PE(shift);
            });
            Vector originEleven = new Vector3D(0,10, 0);
            IMolecule moleculeMOPSix = box.getMoleculeList().get(5);
            moleculeMOPSix.getChildList().forEach(atom -> {
                oldPositionsTwo.get(atom.getIndex()).E(atom.getPosition());
                atom.getPosition().PE(originEleven);
                Vector shift = box.getBoundary().centralImage(atom.getPosition());
                atom.getPosition().PE(shift);
            });
            IMolecule moleculeMOPSeven = box.getMoleculeList().get(6);
            Vector originTen = new Vector3D(10,0, 0);
            moleculeMOPSeven.getChildList().forEach(atom -> {
                oldPositionsTwo.get(atom.getIndex()).E(atom.getPosition());
                atom.getPosition().PE(originTen);
                Vector shift = box.getBoundary().centralImage(atom.getPosition());
                atom.getPosition().PE(shift);
            });
            Vector originEight = new Vector3D(-10,0, 0);
            IMolecule moleculeMOPTwelve = box.getMoleculeList().get(7);
            moleculeMOPTwelve.getChildList().forEach(atom -> {
                oldPositionsTwo.get(atom.getIndex()).E(atom.getPosition());
                atom.getPosition().PE(originEight);
                Vector shift = box.getBoundary().centralImage(atom.getPosition());
                atom.getPosition().PE(shift);
            });*/
        }
      /*  if(ifMOPPresent){
            massMembrane += speciesMOP.getMass();
            box.addNewMolecule(speciesMOP);
            System.out.println(speciesMOP.getMass());
        }*/
     /*   box.setNMolecules(speciesGas, 1);

        Vector originEight = new Vector3D(-8,-8,-8);
        IMolecule moleculeMOPTwelve = box.getMoleculeList().get(1);
        while (oldPositions.size() < moleculeMOPTwelve.getChildList().size()) {
            oldPositions.add(space.makeVector());
        }
        moleculeMOPTwelve.getChildList().forEach(atom -> {
            oldPositions.get(atom.getIndex()).E(atom.getPosition());
            atom.getPosition().PE(originEight);
            Vector shift = box.getBoundary().centralImage(atom.getPosition());
            atom.getPosition().PE(shift);
        });
        System.out.println(moleculeMOPTwelve.getChildList().get(0).getPosition());
*/

        if(ifSecondGasPresent){
            box.addNewMolecule(speciesSecondGas);
        }

        if(ifGraphenePresent && ifMOPPresent){
            if(ifSecondGasPresent){
                sm = new SpeciesManager.Builder().addSpecies(speciesGrapheneOne).addSpecies(speciesMOP).addSpecies(speciesGas).addSpecies(speciesSecondGas).build();
            } else {
                sm = new SpeciesManager.Builder().addSpecies(speciesGrapheneOne).addSpecies(speciesMOP).addSpecies(speciesGas).build();
            }
        }else if (ifGraphenePresent){
            if(ifSecondGasPresent){
                sm = new SpeciesManager.Builder().addSpecies(speciesGrapheneOne).addSpecies(speciesGas).addSpecies(speciesSecondGas).build();
            } else {
                sm = new SpeciesManager.Builder().addSpecies(speciesGrapheneOne).addSpecies(speciesGas).build();
            }
        } else if(ifMOPPresent){
            if(ifSecondGasPresent){
                sm = new SpeciesManager.Builder().addSpecies(speciesMOP).addSpecies(speciesGas).addSpecies(speciesSecondGas).build();
            } else {
                sm = new SpeciesManager.Builder().addSpecies(speciesMOP).addSpecies(speciesGas).build();
            }
        } else {
            if (ifSecondGasPresent) {
                sm = new SpeciesManager.Builder().addSpecies(speciesGas).addSpecies(speciesSecondGas).build();
            } else {
                sm = new SpeciesManager.Builder().addSpecies(speciesGas).build();
            }
        }

       // System.out.println(speciesMOP.getMass());
       // System.exit(1);
       /* ArrayList<ArrayList<Integer>> connectedAtoms1 = pdbReaderMOP.getConnectivity();
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
        Map<String[], List<int[]>> bondTypesMap2 = null;
        Map<String[], List<int[]>> angleTypesMap2 = null;
        Map<String[], List<int[]>> torsionTypesMap2 = null;
        int[] harmonicBond = null;
        List<int[]> bonds = new ArrayList<>();
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

       /* System.out.println(speciesMOP.getMass());
        System.out.println(Gram.UNIT.fromSim(speciesMOP.getMass()));
        Unit cm = new PrefixedUnit(Prefix.CENTI, Meter.UNIT);
        Unit cm3 = new CompoundUnit(new Unit[]{cm}, new double[]{3});
        UnitRatio gcm3 = new UnitRatio(Gram.UNIT, cm3);
        System.out.println(gcm3.toSim(1.082) +" exp");
        double valBox = boxSize.getX(0);
        double mopDen = speciesMOP.getMass()/(valBox*valBox*valBox);
        System.out.println(mopDen);
        System.out.println(gcm3.fromSim(mopDen));*/
   //     System.out.println(gcm3.toSim(1.082)/gcm3.toSim(mopDen));
        //SetBonding Potential
        UniversalSimulation.makeAtomPotentials(sm);
        PotentialMasterBonding pmBonding = new PotentialMasterBonding(sm, box);
        PotentialMasterCell potentialMasterCell = new PotentialMasterCell(getSpeciesManager(), box, 5, pmBonding.getBondingInfo());
        //MOP
        // SetPotential.SetBondStrech(speciesMOP, bondTypesMap1, angleTypesMap1, torsionTypesMap1,bondsNum1,bondList1, quadrupletsSorted1, atomIdentifierMapModified1,atomicPotMap1, pmBonding);
        //GasOne
        if (!isGasTraPPE) {
            SetPotential.setBondStretch(speciesGas, bondTypesMap2, angleTypesMap2, torsionTypesMap2, bondsNum2, bondList2, quadrupletsSorted2, atomIdentifierMapModified2, atomicPotMap2, pmBonding);
        } else {
            String[] parts = confNameGasOne.split("//");
            String gasName = parts[parts.length-1];
            SetPotential.setBondStretchTraPPE(speciesGas, bonds, pmBonding, gasName);
        }

        //GasTwo
        if(ifSecondGasPresent){
            SetPotential.setBondStretch(speciesSecondGas,bondTypesMapSecondGas, angleTypesMapSecondGas, torsionTypesMapSecondGas, bondsNumSecondGas, bondListSecondGas,quadrupletsSortedSecondGas, atomIdentifierMapModifiedSecondGas, atomicPotMapSecondGas, pmBonding);
        }
        List<AtomType> listGas;
        List<List<AtomType>> listMOPGasPairs = null;

        //NonBonded
        potentialMasterCell.doAllTruncationCorrection = false;
        if (ifSecondGasPresent) {
            listGas = SetPotential.listTwoSpeciesPairs(speciesGas.getUniqueAtomTypes(), speciesSecondGas.getUniqueAtomTypes());
        } else {
            listGas = speciesGas.getUniqueAtomTypes();
        }

        if (ifMOPPresent) {
            listMOPGasPairs = SetPotential.listGrapheneSpecial(speciesMOP.getUniqueAtomTypes(), listGas);
        } else {
            listMOPGasPairs = SetPotential.listFinal(listGas);
        }
        if (ifSecondGasPresent) {
            listMOPGasPairs = SetPotential.listGrapheneSpecial(listGas, speciesSecondGas.getUniqueAtomTypes());
        }

        //MOP-Gas

        // System.out.println(listMOPGas);

        LJUFF[] p2LJMOPGas = new LJUFF[listMOPGasPairs.size()];
        P2Electrostatic[] p2ElectroMOPGas = new P2Electrostatic[listMOPGasPairs.size()];
        IPotential2[] p2mopgas = new IPotential2[listMOPGasPairs.size()];
        //   LJCOMPASS[] p2LJMOPGasCOMPASS = new LJCOMPASS[listMOPGasPairs.size()];
        if (isGasTraPPE) {
            SetPotential.doLJElectrostatic(listMOPGasPairs, potentialMasterCell, p2LJMOPGas, p2ElectroMOPGas, listMOPGasPairs.size(), truncCal, doElectrostatics, true);
        } else {
            SetPotential.doLJElectrostatic(listMOPGasPairs, potentialMasterCell, p2LJMOPGas, p2ElectroMOPGas, p2mopgas, listMOPGasPairs.size(), truncCal, doElectrostatics);
        }

        if (ifMOPPresent) {
            List<List<AtomType>> listGasGas = SetPotential.listFinal(speciesGas.getUniqueAtomTypes());
            LJUFF[] p2LJGas = new LJUFF[listGasGas.size()];
            IPotential2[] p2ljGas = new IPotential2[listGasGas.size()];
            P2Electrostatic[] p2ElectroGas = new P2Electrostatic[listGasGas.size()];
            //   LJCOMPASS[] p2LJGasCOMPASS = new LJCOMPASS[listGasGas.size()];
            if (isGasTraPPE) {
                SetPotential.doLJElectrostatic(listGasGas, potentialMasterCell, p2LJGas, p2ElectroGas, listGasGas.size(), truncatedRadius, false, true);
            } else {
                SetPotential.doLJElectrostatic(listGasGas, potentialMasterCell, p2LJGas, p2ElectroGas, p2ljGas, listGasGas.size(), truncatedRadiusLJ, doElectrostatics);
            }
        }
        //Graphene-Gas
        if (ifGraphenePresent) {
            List<List<AtomType>> listGrapheneGasFinal = SetPotential.listGrapheneSpecial(speciesGrapheneOne.getUniqueAtomTypes(), listGas);
            //  List<List<AtomType>> listGrapheneGasFinal = SetPotential.listFinal(listGrapheneGasUnique);
            LJUFF[] p2LJGrapheneGas = new LJUFF[listGrapheneGasFinal.size()];
            //   SetPotential.doLJ(listGrapheneGasFinal, potentialMasterCell,  p2LJGrapheneGas, listGrapheneGasFinal.size(), truncatedRadiusLJ, isGasTraPPE);
        }

        integrator = new IntegratorMC(potentialMasterCell, random, Kelvin.UNIT.toSim(temperature), box);
        potentialMasterCell.doOneTruncationCorrection = true;
        potentialMasterCell.init();
        if (ifMoveRotateMoves) {
            mcMoveMolecule = new MCMoveMolecule(random, potentialMasterCell, box);
            mcMoveMolecule.setStepSize( sigma);
            ((MCMoveStepTracker) mcMoveMolecule.getTracker()).setTunable(false);
            mcMoveMoleculeRotate = new MCMoveMoleculeRotate(random, potentialMasterCell, box);
            mcMoveMoleculeRotate.setStepSize( sigma);
            ((MCMoveStepTracker) mcMoveMoleculeRotate.getTracker()).setTunable(false);
            ((MoleculeSourceRandomMolecule) mcMoveMolecule.getMoleculeSource()).setSpecies(speciesGas);
            ((MoleculeSourceRandomMolecule) mcMoveMoleculeRotate.getMoleculeSource()).setSpecies(speciesGas);
            integrator.getMoveManager().addMCMove(mcMoveMoleculeRotate);
            integrator.getMoveManager().addMCMove(mcMoveMolecule);
            if (makeAllMove) {
                ((MoleculeSourceRandomMolecule) mcMoveMolecule.getMoleculeSource()).setSpecies(speciesMOP);
                ((MoleculeSourceRandomMolecule) mcMoveMoleculeRotate.getMoleculeSource()).setSpecies(speciesMOP);
            }//else {

            // ((MoleculeSourceRandomMolecule) mcMoveMolecule.getMoleculeSource()).setSpecies(speciesMOP);
            // ((MoleculeSourceRandomMolecule) mcMoveMoleculeRotate.getMoleculeSource()).setSpecies(speciesMOP);

            //  }
              ((MoleculeSourceRandomMolecule) mcMoveMolecule.getMoleculeSource()).setSpecies(speciesGas);
            ((MoleculeSourceRandomMolecule) mcMoveMoleculeRotate.getMoleculeSource()).setSpecies(speciesGas);
           //    }
        }
            BoxInflate inflater = new BoxInflate(box, space);
            inflater.actionPerformed();

            MCMoveInsertDelete mcMoveID = new MCMoveInsertDelete(potentialMasterCell, random, space);
            mcMoveID.setBox(box);
            mcMoveID.setMu(mu1);
            mcMoveID.setSpecies(speciesGas);
            integrator.getMoveManager().addMCMove(mcMoveID);

            if (ifSecondGasPresent) {
                MCMoveInsertDelete mcMoveIDTwo = new MCMoveInsertDelete(potentialMasterCell, random, space);
                mcMoveID.setBox(box);
                mcMoveID.setMu(mu1);
                mcMoveID.setSpecies(speciesSecondGas);
                integrator.getMoveManager().addMCMove(mcMoveID);
            }


            //potential = new P2LennardJones(sigma, 1.0);
            //
           /* if (truncCal > 0.5 * box.getBoundary().getBoxSize().getX(0) || truncCal > 0.5 * box.getBoundary().getBoxSize().getX(1) || truncCal > 0.5 * box.getBoundary().getBoxSize().getX(2)) {
                double valX = 0.5 * box.getBoundary().getBoxSize().getX(0);
                double valY = 0.5 * box.getBoundary().getBoxSize().getX(1);
                double valZ =0.5 * box.getBoundary().getBoxSize().getX(2);
                throw new RuntimeException("Truncation radius too large.  Max allowed is " + 0.5 * box.getBoundary().getBoxSize().getX(0) + " " + valX + " "+ valY + " " + valZ +" "  +truncCal);

            }
        double valX = 0.5 * box.getBoundary().getBoxSize().getX(0);
        double valY = 0.5 * box.getBoundary().getBoxSize().getX(1);
        double valZ =0.5 * box.getBoundary().getBoxSize().getX(2);
        System.out.println(0.5 * box.getBoundary().getBoxSize().getX(0) + " " + valX + " "+ valY + " " + valZ +" "  +truncCal);
            System.exit(1);*/

    }

    public static void main(String[] args) throws IOException {
        GCMCMOPParams params = new GCMCMOPParams();
        ParseArgs.doParseArgs(params, args);
        int numAtomOne = params.numAtomOne;
        int numAtomTwo = params.numAtomTwo;
        boolean ifautoMOP = params.ifautoMOP;
        String autoMOPGeom = params.autoMOPGeom;
      //  String confNameOne = params.confNameOne;
       // String confNameGasOne = params.confNameGasOne;
      //  String confNameGasTwo = params.confNameGasTwo;
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
        boolean ifMOPPresent = params.ifMOPPresent;
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
       // System.out.println(mu1 + " mu1");
        List<Integer> muValues = new ArrayList<>();
        List<Integer> temperatureList = new ArrayList<>();

        List<String> mopNames = new ArrayList<>();
        List<String> gasNames = new ArrayList<>();


        List<Double> mopTemp = new ArrayList<>();

       mopNames.add(params.confNameEight);


       // gasNames.add(params.confNamegasThree);
        gasNames.add(params.confNameGasTwo);
        gasNames.add(params.confNamegas);
        gasNames.add(params.confNamegasOne);

      //  temperatureList.add(77);
       // temperatureList.add(77);
        double t1Start = System.nanoTime();
       for (int i = (int) mu1; i>params.muLimit; i= (int) (i+muDecrease)){
            muValues.add(i);
        }
        System.out.println(muValues);
        System.out.println(boxSize);
        //int numSteps = params.numSteps;
     for(int a =0; a<gasNames.size(); a++){
         String confNamegas = gasNames.get(a);
         System.out.println(confNamegas);
      //   double tempActual =temperatureList.get(a);
     //    System.out.println(tempActual);
         for (int p = 0; p < mopNames.size(); p++) {
             // double temperature = mopNames.get(p);
             String confNameMOPOne = mopNames.get(p);
             System.out.println(confNameMOPOne);
             for(int i=0; i<muValues.size(); i++) {
                 int numSteps = params.numSteps;
                 int mu = muValues.get(i);
              //   System.out.println(mu);
                 GCMCMOP sim = new GCMCMOP(confNameMOPOne, confNamegas, confNamegas, confNameGraphene, centreMOP, centreMOPTwo, centreMOPThree, centreMOPFour, grapheneOne, grapheneTwo, grapheneThree, grapheneFour, grapheneFive, grapheneSix, grapheneSeven, grapheneEight, grapheneNine, grapheneTen, grapheneEleven, grapheneTwelve, grapheneThirteen, numAtomOne, numAtomTwo, temperature, truncatedRadius, truncatedRadiusLJ, sigma, mu, mu2, ifGraphenePresent, ifSecondGasPresent, ifMultipleGraphenePresent, ifMoveRotateMoves, boxSize, makeAllMove, doElectrostatics, multiplier, isGasCOMPASS, isGasTraPPE, ifMOPPresent, ifCOFPresent, autoMOPGeom, ifautoMOP);
                 if (massMembrane < 1) {
                     massMembrane = 1;
                 }
                 // System.out.println(Constants.BOLTZMANN_K);

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


                 if (doGraphics) {
                     final SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, "GCMC");
                     DiameterHashByType dhbt = (DiameterHashByType) simGraphic.getDisplayBox(sim.box).getDiameterHash();
                     if (ifMOPPresent) {
                       //  dhbt.setDiameter(sim.speciesMOP.getAtomType(0), 1);
                  /*  dhbt.setDiameter(sim.speciesMOP.getAtomType(1), 1);
                    dhbt.setDiameter(sim.speciesMOP.getAtomType(2), 1);
                    dhbt.setDiameter(sim.speciesMOP.getAtomType(3), 1);
                    dhbt.setDiameter(sim.speciesMOP.getAtomType(4), 1);
                    dhbt.setDiameter(sim.speciesMOP.getAtomType(4), 0.2);*/
                     }

                     if (ifGraphenePresent) {
                         for (int j = 0; i < sim.speciesGrapheneOne.getUniqueAtomTypes().size(); i++) {
                             dhbt.setDiameter(sim.speciesGrapheneOne.getAtomType(j), 1.5);
                         }
                         dhbt.setDiameter(sim.speciesGas.getAtomType(0), 1.5);
                         //dhbt.setDiameter(sim.speciesGas.getAtomType(1), 1);
                         //dhbt.setDiameter(sim.speciesGrapheneOne.getAtomType(0), 1.2);
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
                         ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGrapheneOne.getTypeByName("HK"), ColorExtra.Violet);
                         ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGrapheneOne.getTypeByName("OJ"), ColorExtra.gold);
                         ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGrapheneOne.getTypeByName("OK"), ColorExtra.cornflowerblue);

                     }

                     if (ifGraphenePresent) {
                         ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGrapheneOne.getAtomType(0), Color.gray);
                         ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGrapheneOne.getAtomType(1), Color.gray);
                         ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGrapheneOne.getAtomType(2), ColorExtra.lightcyan);
                         ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGrapheneOne.getAtomType(3), Color.red);
                         ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGrapheneOne.getAtomType(4), Color.gray);
                         ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGrapheneOne.getAtomType(5), ColorExtra.lightcyan);
                         ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGrapheneOne.getAtomType(6), ColorExtra.lightcyan);
                         ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGrapheneOne.getAtomType(7), ColorExtra.lightcyan);
                     }


                     //MOP
                     if (ifMOPPresent) {
                         dhbt.setDiameter(sim.speciesMOP.getAtomType(0), 3);
                  /*  dhbt.setDiameter(sim.speciesMOP.getAtomType(1), 2);
                    dhbt.setDiameter(sim.speciesMOP.getAtomType(2), 1);
                    dhbt.setDiameter(sim.speciesMOP.getAtomType(3), 2);
                    dhbt.setDiameter(sim.speciesMOP.getAtomType(4), 2);
                    dhbt.setDiameter(sim.speciesGas.getAtomType(0), 1.5);*/
                         // dhbt.setDiameter(sim.speciesGas.getAtomType(1), 1);
                         ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesMOP.getAtomType(0), ColorExtra.copper);
                  /*  ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesMOP.getAtomType(1), ColorExtra.red);
                    ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesMOP.getAtomType(2), ColorExtra.blue);
                    ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesMOP.getAtomType(3), ColorExtra.gray);
                    ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesMOP.getAtomType(4), ColorExtra.firebrick);*/
                     }

                     if (ifSecondGasPresent) dhbt.setDiameter(sim.speciesSecondGas.getAtomType(0), 1.5);
                     //if(ifSecondGasPresent)dhbt.setDiameter(sim.speciesSecondGas.getAtomType(1), 1);

                     dhbt.setDiameter(sim.speciesGas.getAtomType(0), 3.5);
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
                 int bs = params.numSteps / (pInterval * 5);
                 if (bs == 0) bs = 1;
             /* MeterPressure pMeter = new MeterPressure(sim.box, sim.integrator.getPotentialCompute());
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
                 int numsBins = params.numBins;
                 DoubleRange numBox = new DoubleRange(-(boxSize.getX(0))/2, (boxSize.getX(0))/2);
                 if(params.doHistogram){
                     System.out.println("Numbins " + numsBins);
                     HistogramVectorSimple histogramVectorSimple = new HistogramVectorSimple(numsBins,numsBins,numsBins, numBox, numBox, numBox, sim.speciesGas, sim.box, numSteps);
                     histogramVectorSimple.actionPerformed();
                 }

                 //     long[][][] count = histogramVectorSimple.getCounts();
                 // System.out.println(Arrays.deepToString(count));

                 AccumulatorAverage densityAccumulator = new AccumulatorAverageFixed(bs);
                 densityMeter.setSpecies(sim.speciesGas);
                 DataPumpListener pumpDensity = new DataPumpListener(densityMeter, densityAccumulator);
                 sim.integrator.getEventManager().addListener(pumpDensity);

                 if (ifSecondGasPresent) {
                     MeterDensity densityMeterGasTwo = new MeterDensity(sim.box);
                     AccumulatorAverage densityAccumulatorGastwo = new AccumulatorAverageFixed(bs);
                     densityMeterGasTwo.setSpecies(sim.speciesSecondGas);
                     DataPumpListener pumpDensityGasTwo = new DataPumpListener(densityMeterGasTwo, densityAccumulatorGastwo);
                     sim.integrator.getEventManager().addListener(pumpDensityGasTwo);
                 }
                 if(params.doHistogram){
                     final Vector[] vec = {new Vector3D()};
                     IMoleculeList molecules = sim.box.getMoleculeList(sim.speciesGas);
                     final HistogramVectorSimple targHist = new HistogramVectorSimple(numsBins,numsBins,numsBins,numBox, numBox, numBox, sim.speciesGas, sim.box , numSteps);

                     IntegratorListener histListenerTarget = new IntegratorListener() {
                         public void integratorStepStarted(IntegratorEvent e) {}

                         public void integratorStepFinished(IntegratorEvent e) {
                             for(int v=0; v<molecules.size(); v++){
                                 vec[0] = molecules.get(v).getChildList().get(1).getPosition();
                                 targHist.addValue(vec[0].getX(0),vec[0].getX(1), vec[0].getX(2));
                             }

                         }

                         public void integratorInitialized(IntegratorEvent e) {}
                     };


                     System.out.println("collecting histograms");


                     // only collect the histogram if we're forcing it to run the reference system
                     sim.integrator.getEventManager().addListener(histListenerTarget);
                 }
                 long t1 = System.currentTimeMillis();
                 sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, params.numSteps));
                 long t2 = System.currentTimeMillis();
                 //System.out.println("mu1 " + mu1);
                 System.out.println("runtime: " + (t2 - t1) * 0.001);
          /*     Unit pUnit = Bar.UNIT;
           //     Unit MegaPascal = new PrefixedUnit(Prefix.MEGA, Pascal.UNIT);
               Unit kiloPascal = new PrefixedUnit(Prefix.KILO, Pascal.UNIT);
                double avgP = pAccumulator.getData(pAccumulator.AVERAGE).getValue(0);
                double errP = pAccumulator.getData(pAccumulator.ERROR).getValue(0);
                double corP = pAccumulator.getData(pAccumulator.BLOCK_CORRELATION).getValue(0);
             //   System.out.println("P (MPa) " + MegaPascal.fromSim(avgP) + " MPa " + MegaPascal.fromSim(errP) + " " + corP);
                 System.out.println("P (kPa) " + kiloPascal.fromSim(avgP) + " kPa " + kiloPascal.fromSim(errP) + " " + corP);
                System.out.println("P (Bar) " + pUnit.fromSim(avgP) + " bar " + pUnit.fromSim(errP) + " " + corP);*/
                 //double muR = Constants.BOLTZMANN_K * temperature * Math.log(avgP/temperature);
                 //System.out.println("muR: " + muR);

                 double avgRho = densityAccumulator.getData(densityAccumulator.AVERAGE).getValue(0);
                 double errRho = densityAccumulator.getData(densityAccumulator.ERROR).getValue(0);
                 double corRho = densityAccumulator.getData(densityAccumulator.BLOCK_CORRELATION).getValue(0);
                 Unit dm = new PrefixedUnit(Prefix.DECI, Meter.UNIT);
                 Unit dm3 = new CompoundUnit(new Unit[]{dm}, new double[]{3});
                 Unit m3 = new CompoundUnit(new Unit[]{Meter.UNIT}, new double[]{3});
                 Unit kgm3 = new UnitRatio(new PrefixedUnit(Prefix.KILO, Gram.UNIT), m3);
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
                 //  System.out.println("Actual "+ avgRho);
                 //System.out.println("nm " + nm3.fromSim(sim.box.getBoundary().volume()));
                System.out.println("Num atoms : " + numAtomsAvg + " " + errnum + " " + cornum + " " + (sim.speciesGas.getMass()*numAtomsAvg*1000/sim.speciesMOP.getMass()));
                 // System.out.println("num atoms (/cm3) "+cm3.fromSim( sim.box.getBoundary().volume())*errRho + " " + cm3.fromSim( sim.box.getBoundary().volume())*errRho);
                 // System.out.println("num atoms (/m3) " +m3.fromSim( sim.box.getBoundary().volume())*errRho + " " + m3.fromSim( sim.box.getBoundary().volume())*errRho);
                 // System.out.println("num atoms (/nm3) " +nm3.fromSim( sim.box.getBoundary().volume())*errRho + " " + nm3.fromSim( sim.box.getBoundary().volume())*errRho);
                 // System.out.println("num atoms (/nm3) " +nm_3.fromSim(avgRho) + " " + nm_3.fromSim(errRho) );
                 // System.out.println("Rho (mol/dm3) "+kgm3.fromSim(avgRho)/sim.speciesGas.getMass() + " " + sim.speciesGas.getMass());
                 //AtomType C1 = new AtomType(Carbon.INSTANCE, "C1");
//            System.out.println("Amount sorbed " + numAtomsAvg*1000*sim.speciesGas.getMass() /(sim.speciesMOP.getMass()) + " err " +errnum*1000*sim.speciesGas.getMass() /(sim.speciesMOP.getMass()));


                 massGas = numAtomsAvg * sim.speciesGas.getMass();
                 // System.out.println(massGas + " Graphene " + massGraphene  + " single graphene " + sim.speciesGrapheneOne.getMass());
                 //  System.out.println("Mass ratio " + massGas/(massGraphene));
                 //   System.out.println(massMOP + " "+massGraphene );
                 //  System.out.println("Amount sorbed " + 1000*massGas/(massGraphene+massMOP)  + " err " + 1000*errnum*sim.speciesGas.getMass()/(massGraphene+massMOP));
                 //      System.out.println("Amount sorbed per Graphene basis " + massGas/((massGraphene))  + " err " +errnum*sim.speciesGas.getMass()/(massGraphene));
                 //  System.out.println("Density (g/cm3) " + gcm3.fromSim(avgRho));
//          System.out.println("Mass fraction : "+2*sim.speciesMOP.getMass()*100/(sim.speciesGrapheneOne.getMass()*12));
          /*  double avgPE = energyAccumulator.getData(energyAccumulator.AVERAGE).getValue(0) / numAtomsAvg;
            double errPE = energyAccumulator.getData(energyAccumulator.ERROR).getValue(0) / numAtomsAvg;
            double corPE = energyAccumulator.getData(energyAccumulator.BLOCK_CORRELATION).getValue(0);*/
                 //  System.out.println("PE " + avgPE + " " + errPE + " " + corPE);
                 // double numAtomGasOne = nMoleculesGasOne.getDataAsScalar();
                 //System.out.println("Atoms gas One: " +numAtomGasOne);

                 if (ifSecondGasPresent) {
                     double numAtomGasTwo = nMoleculesGasTwo.getDataAsScalar();
                     System.out.println("Atoms gas Two: " + numAtomGasTwo);
                 }
//        System.out.println(numAtomsAvg*sim.speciesGas.getMass()*100/(sim.speciesGrapheneOne.getMass()) + " " +numAtomsAvg*sim.speciesGas.getMass()+ " "+sim.speciesGrapheneOne.getMass());
            //     System.out.println("rho (kg/m3) " + kgm3.fromSim(avgRho) + " " + kgm3.fromSim(errRho) + " " + corRho);
                 System.out.println("rho " + (avgRho) + " " + errRho + " " + corRho + " " + (errRho / avgRho) * 100);
                 // System.out.println("uptake  " + 22.414 * kgm3.fromSim(avgRho) / sim.speciesGas.getMass());
                 // System.out.println("Rho (mol/dm3) " + moldm3.fromSim(avgRho) + " " +errRho + " " +corRho);
                 // expected values based on 10^8 steps
                 // stdev based on 10^8 steps for N=500 uncertainty, scaled up to 10^6 steps
                 //   other sims have significant correlation
                 // 4 sigma should fail 1 in 16,000 runs

                 double expectedP = 0.0815; // finite size effect smaller than uncertainty
                 double stdevP = 0.03;
       /* if (Double.isNaN(avgP) || Math.abs(avgP - expectedP) / stdevP > 4) {
            System.exit(1);
        }*/


                 double expectedPE = -4.492; // finite size effect comparable to uncertainty
                 double stdevPE = 0.09;
                 if(params.doHistogram){
                     double[][][] count = targHist.getHistogram();
                     System.out.println(Arrays.deepToString(count));
                 }



           /* if (Double.isNaN(avgPE) || Math.abs(avgPE - expectedPE) / stdevPE > 4) {
                System.exit(2);
            }
            double[] boxDimension = boxSize.toArray();
            // System.out.println(sim.speciesGrapheneOne.getMass()*2/(boxDimension[0]* boxDimension[1]*boxDimension[2]));
            double expectedRho = 0.65; // finite size effect smaller than uncertainty
            double stdevRho = 0.006;

            if (Double.isNaN(avgRho) || Math.abs(avgRho - expectedRho) / stdevRho > 4) {
                System.exit(3);
            }*/
                 //  System.out.println("\n");
             }
             System.out.println("\n");
         }


         System.out.println("\n");
     }

        double t2 = System.nanoTime();
        System.out.println((t2-t1Start)/Math.pow(10,9) );


    }

    public void getPolyhedra(Map<Double, List<Integer[]>> distMap,List<Double> slopeXY, List<Double> angleZ, String geometry){
        List<Vector> listPositionsUnitPoly = new ArrayList<>();
        switch (geometry){
            case "tetra":
                listPositionsUnitPoly = getTetrahedron();
                break;
            case "rhombic":
                listPositionsUnitPoly = getRhombic();
                break;
            case "cubic":
                listPositionsUnitPoly = getCubic();
                break;
            case "icosa":
                listPositionsUnitPoly = getIcosahedron();
                break;
            case "dodeca":
                listPositionsUnitPoly =getDodecahedron();
                break;
        }
        edgesPolyhedra(distMap, listPositionsUnitPoly);
        printResults(distMap, listPositionsUnitPoly,  slopeXY, angleZ);
    }

    public List<Vector> getTetrahedron(){
        List<Vector> listVectorTetra = new ArrayList<>();
        Vector vec1 = Vector.of(1, 1, 1);
        Vector vec2 = Vector.of(-1, -1, 1);
        Vector vec3 = Vector.of(-1, 1, -1);
        Vector vec4 = Vector.of(1, -1, -1);
        listVectorTetra.add(vec1);
        listVectorTetra.add(vec2);
        listVectorTetra.add(vec3);
        listVectorTetra.add(vec4);
        return listVectorTetra;
    }

    public List<Vector> getRhombic(){
        List<Vector> listVecRhombic = new ArrayList<>();
        Vector vec1 = Vector.of(1,0 , 0);
        Vector vec2 = Vector.of(-1,0 ,0 );
        Vector vec3 = Vector.of(0,1 , 0);
        Vector vec4 = Vector.of(0, -1, 0);
        Vector vec5 = Vector.of(0,0 , 1);
        Vector vec6 = Vector.of(0,0 , -1);
        listVecRhombic.add(vec1);
        listVecRhombic.add(vec2);
        listVecRhombic.add(vec3);
        listVecRhombic.add(vec4);
        listVecRhombic.add(vec5);
        listVecRhombic.add(vec6);
        return listVecRhombic;
    }

    public List<Vector> getCubic(){
        List<Vector> listVecCubic = new ArrayList<>();
        Vector vec1 = Vector.of(1, 1, 1);
        Vector vec2 = Vector.of(-1, -1, -1);
        Vector vec3 = Vector.of(-1, 1, 1);
        Vector vec4 = Vector.of(1, -1, 1);
        Vector vec5 = Vector.of(1, 1, -1);
        Vector vec6 = Vector.of(-1, -1, 1);
        Vector vec7 = Vector.of(1, -1, -1);
        Vector vec8 = Vector.of(-1, 1, -1);
        listVecCubic.add(vec1);
        listVecCubic.add(vec2);
        listVecCubic.add(vec3);
        listVecCubic.add(vec4);
        listVecCubic.add(vec5);
        listVecCubic.add(vec6);
        listVecCubic.add(vec7);
        listVecCubic.add(vec8);
        return listVecCubic;
    }

    public List<Vector> getIcosahedron(){
        double phi = (1+ Math.sqrt(5))/2;
        double phi2 = phi*phi;
        List<Vector> listVecIcosa = new ArrayList<>();
        Vector vec1 = Vector.of(0,1 , phi);
        Vector vec2 = Vector.of(0, 1, -phi);
        Vector vec3 = Vector.of(0, -1, phi);
        Vector vec4 = Vector.of(0, -1, -phi);
        Vector vec5 = Vector.of(1,phi ,0 );
        Vector vec6 = Vector.of(1,-phi,0 );
        Vector vec7 = Vector.of(-1, phi,0 );
        Vector vec8 = Vector.of(-1,-phi,0 );
        Vector vec9 = Vector.of(phi, 0,1 );
        Vector vec10 = Vector.of(-phi,0 ,1 );
        Vector vec11 = Vector.of(phi,0 ,-1 );
        Vector vec12 = Vector.of(-phi, 0, -1);
        listVecIcosa.add(vec1);
        listVecIcosa.add(vec2);
        listVecIcosa.add(vec3);
        listVecIcosa.add(vec4);
        listVecIcosa.add(vec5);
        listVecIcosa.add(vec6);
        listVecIcosa.add(vec7);
        listVecIcosa.add(vec8);
        listVecIcosa.add(vec9);
        listVecIcosa.add(vec10);
        listVecIcosa.add(vec11);
        listVecIcosa.add(vec12);
        return listVecIcosa;
    }

    public List<Vector> getDodecahedron(){
        double phi = (1+ Math.sqrt(5))/2;
        double phi2 = phi*phi;
        List<Vector> listVecDodeca = new ArrayList<>();
        Vector vec1 = Vector.of(1, 0, phi2);
        Vector vec2 = Vector.of(1, 0, -phi2);
        Vector vec3 = Vector.of(-1, 0, phi2);
        Vector vec4 = Vector.of(-1, 0, -phi2);
        Vector vec5 = Vector.of(phi, phi,phi );
        Vector vec6 = Vector.of(-phi, -phi,-phi );
        Vector vec7 = Vector.of(phi, -phi, -phi);
        Vector vec8 = Vector.of(phi, phi, -phi);
        Vector vec9 = Vector.of(-phi, phi,-phi );
        Vector vec10 = Vector.of(-phi,phi ,phi );
        Vector vec11 = Vector.of(-phi,-phi , phi);
        Vector vec12 = Vector.of(phi,-phi ,phi);
        Vector vec13 = Vector.of(0, phi2,1 );
        Vector vec14 = Vector.of(0, phi2, -1);
        Vector vec15 = Vector.of(0,-phi2 ,1 );
        Vector vec16 = Vector.of(0, -phi2, -1);
        Vector vec17 = Vector.of(phi2, 1, 0);
        Vector vec18 = Vector.of(phi2, -1, 0);
        Vector vec19 = Vector.of(-phi2, 1, 0);
        Vector vec20 = Vector.of(-phi2, -1, 0);
        listVecDodeca.add(vec1);
        listVecDodeca.add(vec2);
        listVecDodeca.add(vec3);
        listVecDodeca.add(vec4);
        listVecDodeca.add(vec5);
        listVecDodeca.add(vec6);
        listVecDodeca.add(vec7);
        listVecDodeca.add(vec8);
        listVecDodeca.add(vec9);
        listVecDodeca.add(vec10);
        listVecDodeca.add(vec11);
        listVecDodeca.add(vec12);
        listVecDodeca.add(vec13);
        listVecDodeca.add(vec14);
        listVecDodeca.add(vec15);
        listVecDodeca.add(vec16);
        listVecDodeca.add(vec17);
        listVecDodeca.add(vec18);
        listVecDodeca.add(vec19);
        listVecDodeca.add(vec20);
        return listVecDodeca;
    }

    public Map<Double, List<Integer[]>> edgesPolyhedra ( Map<Double, List<Integer[]>> distMap, List<Vector> listVect){
        List<Integer[]> doubleList = new ArrayList<>();
        Integer[] points = new Integer[2];
        for (int i = 0; i < listVect.size(); i++) {
            Vector pointA = listVect.get(i);
            for (int j = i + 1; j < listVect.size(); j++) { // Start from i + 1 to avoid duplicates
                Vector pointB = listVect.get(j);

                double distance = calculateDistance(pointA, pointB);
                doubleList = distMap.get(distance);
                if(doubleList == null){
                    doubleList = new ArrayList<>();
                    points = new Integer[]{i, j};
                    doubleList.add(points);
                    distMap.put(distance, doubleList);
                }else {
                    points = new Integer[]{i, j};
                    doubleList.add(points);
                    distMap.put(distance, doubleList);
                }
                //   System.out.println("Distance between point " + i + " and point " + j + ": " + distance);
            }
        }
        return distMap;
    }
    public double calculateDistance(Vector pointA, Vector pointB){
        return pointA.Mv1Squared(pointB);
    }

    public void printResults(Map<Double, List<Integer[]>> distMap,  List<Vector> listVecRhombic,  List<Double> slopeXY, List<Double> angleZ){
        Double smallestKey = null;
        System.out.println(  distMap.keySet() );
        //  System.exit(1);
        double smallestSize = Integer.MAX_VALUE;

        for (Map.Entry<Double, List<Integer[]>> entry : distMap.entrySet()) {
            double currentSize = entry.getKey();
            // Check for smaller size or the same size with a smaller key
            if (currentSize < smallestSize || (currentSize == smallestSize && (smallestKey == null || entry.getKey() < smallestKey))) {
                smallestSize = currentSize;
                smallestKey = entry.getKey();
            }
        }

        if (smallestKey != null) {
            System.out.println("Key with the smallest list size: " + smallestKey);
        } else {
            System.out.println("Map is empty.");
        }
        // System.out.println(distMap.entrySet());
        List<Integer[]> valueList = distMap.get(smallestKey);
        printVal(valueList, listVecRhombic, slopeXY, angleZ);

       /* for (Map.Entry<Double, List<Integer[]>> entry : distMap.entrySet()) {
            Double key = entry.getKey();
            List<Integer[]> valueList = entry.getValue();
            if (valueList.size() > 28) {
            //    System.out.println(key + " : " + Arrays.deepToString(valueList.toArray()) + " Size " + valueList.size());
                printVal(valueList, listVecRhombic, slopeXY, angleZ);
            }
        }*/
    }
    public void printVal( List<Integer[]> valueList,  List<Vector> listVecRhombic, List<Double> slopeXY, List<Double> angleZ){
        for(int i=0; i<valueList.size(); i++){
            System.out.println(Arrays.toString(valueList.get(i)) + " " + valueList.get(i)[0] + " "+ valueList.get(i)[1]   +" " + listVecRhombic.get(valueList.get(i)[0]) + " "+ listVecRhombic.get(valueList.get(i)[1]));
            double slope = slopeCalc(listVecRhombic.get(valueList.get(i)[0]), listVecRhombic.get(valueList.get(i)[1]));
            double zAngle = angleZCalc(listVecRhombic.get(valueList.get(i)[0]), listVecRhombic.get(valueList.get(i)[1]));
            slopeXY.add(slope);
            angleZ.add(zAngle);
        }
        // System.exit(1);
    }
    public double angleZCalc(Vector v1, Vector v2){
        double slope =0;
        double x1= v1.getX(0), x2 = v2.getX(0), y1 = v1.getX(1), y2 = v2.getX(1), z1 = v1.getX(2), z2 = v2.getX(2);
        Vector xNew = new Vector3D(x1, y1, z2);
        Vector y = new Vector3D(x2, y2, z2);
        Vector x = new Vector3D(x1, y1, z1);
        Vector xSub = new Vector3D();
        Vector ySub = new Vector3D();
        xSub.Ev1Mv2(xNew, x);
        ySub.Ev1Mv2(x, y);
        double dotVal = xSub.dot(ySub);
        double mag = Math.sqrt(xSub.squared() * ySub.squared());
        slope = Math.acos(dotVal/mag);
        Degree.UNIT.fromSim(slope);
        System.out.println(Degree.UNIT.fromSim(slope));
        return slope;
    }
   /* public double slp(double x1, double y1, double x2, double y2){
        return (y2-y1)/(x2-x1);
    }*/

    public double  slopeCalc(Vector v1, Vector v2){
        double slope =0;
        double x1= v1.getX(0), x2 = v2.getX(0), y1 = v1.getX(1), y2 = v2.getX(1);
        slope = (y2-y1)/(x2-x1);
        return slope;
    }
    // new Coordinates in XY plane
    public List<Vector> rotatePoint(Vector vectorA, IMolecule molecule, List<Double> theta) {
        List<Vector> newPosn = new ArrayList<>();
        for(int i=1; i<molecule.getChildList().size(); i++ ){
            Vector vec = molecule.getChildList().get(i).getPosition();
            double thetaNew = theta.get(i);
            double xB = vec.getX(0) - vectorA.getX(0);
            double yB = vec.getX(1) - vectorA.getX(1);

            double xBrot = xB * Math.cos(thetaNew) - yB * Math.sin(thetaNew);
            double yBrot = xB * Math.sin(thetaNew) + yB * Math.cos(thetaNew);
            newPosn.add(new  Vector3D(vectorA.getX(0) + xBrot, vectorA.getX(1) + yBrot, 0));
        }
        return newPosn;
    }

    // new Z coordinate
    public List<Vector> calculateZ(List<Vector> newPosn, List<Double> thetaZ) {
        List<Vector> finalPosn = new ArrayList<>();
        for(int i=0; i<thetaZ.size(); i++){
            double x = newPosn.get(i).getX(0);
            double y = newPosn.get(i).getX(0);
            finalPosn.add(new Vector3D(x, y, Math.sqrt(x * x + y * y) * Math.tan(thetaZ.get(i))));
        }
       return finalPosn;
    }

    //left most atom in Box
    public int leftSideMost(IMolecule molecule){
        double minX = Double.POSITIVE_INFINITY;
        int iNum = 0;
        for (int i=0; i<molecule.getChildList().size(); i++){
            Vector vec = molecule.getChildList().get(i).getPosition();
            if (vec.getX(0) < minX){
                minX = vec.getX(0) ;
                iNum = i;
            }
        }
        return iNum;
    }

    //right most atom in Box
    public int rightSideMost(IMolecule molecule){
        double maxX = Double.POSITIVE_INFINITY;
        int iNum = 0;
        for (int i=0; i<molecule.getChildList().size(); i++){
            Vector vec = molecule.getChildList().get(i).getPosition();
            if ((vec.getX(0)  > maxX)){
                maxX = vec.getX(0) ;
                iNum = i;
            }
        }
        return iNum;
    }

    //back most atom in Box
    public int backSideMost(IMolecule molecule){
        double minY = Double.POSITIVE_INFINITY;
        int iNum = 0;
        for (int i=0; i<molecule.getChildList().size(); i++){
            Vector vec = molecule.getChildList().get(i).getPosition();
            if (vec.getX(1) < minY){
                minY = vec.getX(1) ;
                iNum = i;
            }
        }
        return iNum;
    }

    //front most atom in Box
    public int frontSideMost(IMolecule molecule){
        double maxY = Double.POSITIVE_INFINITY;
        int iNum = 0;
        for (int i=0; i<molecule.getChildList().size(); i++){
            Vector vec = molecule.getChildList().get(i).getPosition();
            if ((vec.getX(1)  > maxY)){
                maxY = vec.getX(1) ;
                iNum = i;
            }
        }
        return iNum;
    }

    //bottom most atom in Box
    public int bottomSideMost(IMolecule molecule){
        double minZ = Double.POSITIVE_INFINITY;
        int iNum = 0;
        for (int i=0; i<molecule.getChildList().size(); i++){
            Vector vec = molecule.getChildList().get(i).getPosition();
            if (vec.getX(2) < minZ){
                minZ = vec.getX(2) ;
                iNum = i;
            }
        }
        return iNum;
    }
    //top most atom in Box
    public int topSideMost(IMolecule molecule){
        double maxZ = Double.POSITIVE_INFINITY;
        int iNum = 0;
        for (int i=0; i<molecule.getChildList().size(); i++){
            Vector vec = molecule.getChildList().get(i).getPosition();
            if ((vec.getX(2)  > maxZ)){
                maxZ = vec.getX(2) ;
                iNum = i;
            }
        }
        return iNum;
    }
    public Vector boxSizeVector (Map<Integer, Vector> positions){
        double minX = Double.POSITIVE_INFINITY;
        double maxX = Double.NEGATIVE_INFINITY;
        double minY = Double.POSITIVE_INFINITY;
        double maxY = Double.NEGATIVE_INFINITY;
        double minZ = Double.POSITIVE_INFINITY;
        double maxZ = Double.NEGATIVE_INFINITY;
        Vector boxSize = new Vector3D();
        for (Vector vec : positions.values()) {

            if (vec.getX(0) < minX) minX = vec.getX(0) ;
            if (vec.getX(0)  > maxX) maxX = vec.getX(0) ;
            if (vec.getX(1) < minY) minY = vec.getX(1);
            if (vec.getX(1) > maxY) maxY = vec.getX(1);
            if (vec.getX(2) < minZ) minZ = vec.getX(2);
            if (vec.getX(2) > maxZ) maxZ = vec.getX(2);
        }
        //   System.out.println(minX + " " + minY + " " + minZ);
        //  System.out.println(maxX + " " + maxY + " " + maxZ);
        double valX = (maxX - minX) * 1.1;
        double valY = (maxY - minY) * 1.1;
        double valZ = (maxZ - minZ) * 1.1;
        boxSize = Vector.of(valX, valY, valZ);
        return boxSize;
    }

   public Integrator getIntegrator(){
       return integrator;
   }
    public static class GCMCMOPParams extends ParameterBase {
        public Vector grapheneThirteen = new Vector3D(15.0,15.0, 60.0);
        public Vector grapheneOne = new Vector3D(15.0,15.0, 0.0);
        public Vector grapheneTwo = new Vector3D(0,0, -38.5);
        public Vector grapheneSix = new Vector3D(0,0, -19.5);
        public Vector grapheneThree = new Vector3D(0,0, 19.5);
        public Vector grapheneSeven = new Vector3D(0,0, 38.5);
        public Vector grapheneFour = new Vector3D(15.0,15.0, 30.0);
        public Vector grapheneNine = new Vector3D(15.0,15.0, 40.0);
        public Vector grapheneTen = new Vector3D(15.0,15.0, 50.0);
        public Vector grapheneFive = new Vector3D(15.0,15.0, -10.0);
        public Vector grapheneEight = new Vector3D(15.0,15.0, -40.0);
        public Vector grapheneEleven = new Vector3D(15.0,15.0, -50.0);
        public Vector grapheneTwelve = new Vector3D(15.0,15.0, -60.0);
        public double mu2 = -1900;
        public Vector centreMOPTwo = new Vector3D(0.0,0.0,30);
        public Vector centreMOP = new Vector3D(0.0,0.0,-30);
        public Vector centreMOPThree = new Vector3D(0.0,0.0,-50);
        public Vector centreMOPFour = new Vector3D(0.0,0.0,50);

        public int numAtomOne = 1;
        public int numAtomTwo = 1;

        public double sigma =2.7;
        public int numSteps = 10000000;
        public boolean makeAllMove = false;
        public boolean doHistogram = false;
        //public double partial = -82.2479;
        public double temperature = 330;
        public String confNameGasOne = "F://Avagadro//molecule//ethane" ;
        public String confNamegas = "F://Avagadro//molecule//ethane" ;
        public String confNamegasOne = "F://Avagadro//molecule//ethene" ;
        public String confNameGasTwo ="F://Avagadro//molecule//ch4";

        public String confNamegasThree = "F://Avagadro//molecule//ethane" ;

        public double mu1 = -2400;
        public double muDecrease = -50;
        public int muLimit = -2440;
        public boolean ifCOF = false;
        public String confNameCOF = "F://Avagadro//mop//tetra_CuCOOCH3";
        public String confNameGraphene = "F://Avagadro//holeGO60";
        public String confNamemopFive = "F://Avagadro//mop//tetra_cu_8C";
        public String confNamemopSix = "F://Avagadro//MOP_new//icosahedron//propyl_Co_MOP";
        public Vector boxSize = new Vector3D(40, 35, 42);
     //   public Vector boxSize = new Vector3D(16, 16, 16);
        public double truncatedRadiusLJ = boxSize.getX(0)/2;
        public double truncatedRadius = boxSize.getX(0)/2;
        public int numBins = 6;
        public String confNameEight = "F://Avagadro//mop//newMOPs//tetra_cu_(C6H5)2(CO2)2";
        public boolean doGraphics = true;
        public boolean ifautoMOP = true;
        public String autoMOPGeom ="tetra";
        public boolean isGasCOMPASS = false;
        public boolean isGasTraPPE = true;
        public boolean ifGraphenePresent = false;
        public boolean ifMultipleGraphenePresent = false;
        public boolean ifSecondGasPresent = false;
        public boolean ifMOPPresent = true;
        public boolean ifMoveRotateMoves = true;
        public boolean doElectrostatics = false;
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
