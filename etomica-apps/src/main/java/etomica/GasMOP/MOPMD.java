package etomica.GasMOP;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.*;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.data.*;
import etomica.data.history.HistoryCollapsingAverage;
import etomica.data.history.HistoryCollapsingDiscard;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.meter.MeterPressure;
import etomica.data.meter.MeterTemperature;
import etomica.data.types.DataDouble;
import etomica.graphics.*;
import etomica.integrator.*;
import etomica.integrator.mcmove.MCMoveInsertDelete;
import etomica.integrator.mcmove.MCMoveMolecule;
import etomica.integrator.mcmove.MCMoveMoleculeRotate;
import etomica.lattice.LatticeCubicFcc;
import etomica.molecule.IMolecule;
import etomica.molecule.MoleculeSource;
import etomica.molecule.MoleculeSourceRandomMolecule;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.*;
import etomica.potential.UFF.*;
import etomica.potential.compute.PotentialCompute;
import etomica.potential.compute.PotentialComputeAggregate;
import etomica.simulation.Simulation;
import etomica.simulation.prototypes.MeterTorsionAngle;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesManager;
import etomica.units.*;
import etomica.units.dimensions.Null;
import etomica.util.Constants;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.util.collections.IntArrayList;
import etomica.util.random.IRandom;
import etomica.util.random.RandomMersenneTwister;

import java.awt.*;
import java.util.*;
import java.util.List;

public class MOPMD extends Simulation{
    public IntegratorMC integratorMC;
    public IMolecule molecule;
    public MoleculeSource moleculeSource;

    public IntegratorVelocityVerlet integratorMD;
    public IntegratorListenerNHC nhc;
    public MCMoveMolecule mcMoveMolecule;
    public MCMoveMoleculeRotate mcMoveMoleculeRotate;
    public MCMoveInsertDelete mcMoveID;
    public PotentialMaster potentialMaster;
    public PotentialCompute pcAgg, pcAggMC;
    public ISpecies speciesMOP, speciesGasTwo,  speciesGas, speciesGrapheneOne, speciesGrapheneTwo, speciesGrapheneThree, speciesGrapheneFour;
    public Box box;
    public P2LennardJones potential;
    public SpeciesManager sm;
    public int i =0;
    public Vector translationVector, oldPosition;
    public IAtom atom;
    double uOld;
    double uNew = Double.NaN;
    public AtomSource atomSource;
    boolean fixOverlap;
    Space space;

    public MOPMD(String confNameOne, String confNameTwo,String confNameGasTwo,int numGasOne, int numGasTwo, Vector centreMOP,Vector centreMOPTwo, String confNameGraphene, Vector grapheneOne, Vector grapheneTwo, Vector grapheneThree, Vector grapheneFour,int temperature, int truncatedRadiusCharge, int truncatedRadiusLJ, double sigma, double mu, boolean ifGraphenePresent, boolean isDynamic, int numMOP, int numGraphene, int numGas, Vector boxsize, boolean setMOPInfinite, boolean ifGrapheneFixed, boolean ifMOPPresent, boolean ifMOPFixed, boolean ifSecondGasPresent, boolean ifGOMOPMove) {
        super(Space3D.getInstance());
        PDBReaderMOP pdbReaderMOP = new PDBReaderMOP();
        PDBReaderMOP pdbReaderReplicaNew = new PDBReaderMOP();
        PDBReaderMOP pdbReaderReplicaGasNew = new PDBReaderMOP();
        GeneralGrapheneReader grapheneReader = new GeneralGrapheneReader();
        GeneralGrapheneReader grapheneReaderTwo = new GeneralGrapheneReader();
      //  System.out.println(Constants.BOLTZMANN_K);
       // System.out.println(1/(Constants.BOLTZMANN_K*temperature));
        List<ISpecies> listSpecies = new ArrayList<>();
     //   System.out.println(Constants.BOLTZMANN_K);
     //   System.out.println(1/(Constants.BOLTZMANN_K*temperature));
        if(ifMOPPresent){
            speciesMOP = pdbReaderMOP.getSpecies(confNameOne, true, centreMOP, ifMOPFixed, true);
            addSpecies(speciesMOP);
            listSpecies.add(speciesMOP);
        }
        Vector neworiginGas = new Vector3D(5,0,0);
        if(ifGraphenePresent){
            speciesGrapheneOne = grapheneReader.getSpecies(confNameGraphene, grapheneOne, true);
            addSpecies(speciesGrapheneOne);
            listSpecies.add(speciesGrapheneOne);
            speciesGrapheneTwo = grapheneReader.getSpecies(confNameGraphene, grapheneTwo, true);
            addSpecies(speciesGrapheneTwo);
            listSpecies.add(speciesGrapheneTwo);
            //speciesGrapheneThree = grapheneReader.getSpecies(confNameGraphene, grapheneThree, true);
           // addSpecies(speciesGrapheneThree);
           // listSpecies.add(speciesGrapheneThree);
           // speciesGrapheneFour = grapheneReader.getSpecies(confNameGraphene, grapheneFour, true);
            //addSpecies(speciesGrapheneFour);
            //listSpecies.add(speciesGrapheneFour);
        }
        speciesGas = pdbReaderReplicaNew.getSpeciesGas(confNameTwo, true, neworiginGas, false);
        addSpecies(speciesGas);
        listSpecies.add(speciesGas);
        if(ifSecondGasPresent){
            speciesGasTwo = pdbReaderReplicaGasNew.getSpeciesGas(confNameGasTwo, true, neworiginGas, false);
            addSpecies(speciesGasTwo);
            listSpecies.add(speciesGasTwo);
        }
        //speciesMOP
        SetPotential doPotential = new SetPotential();
        List<List<AtomType>> pairAtoms1 = new ArrayList<>();
        List<List<AtomType>> pairAtomsGraphene = new ArrayList<>();
        int pairAtomSize1 =0, pairAtomsGrapheneSize=0;
        if (ifMOPPresent){
            pairAtoms1 = doPotential.getSpeciesPairs(speciesMOP);
            pairAtomSize1 = pairAtoms1.size();
        }
        if(ifGraphenePresent && !ifGrapheneFixed){
            pairAtomsGraphene = doPotential.getSpeciesPairs(speciesGrapheneOne);
            pairAtomsGrapheneSize = pairAtomsGraphene.size();
        }

        //species Gas
        List<List<AtomType>> pairsAtoms2 = doPotential.getSpeciesPairs(speciesGas);
        int pairAtomSize2 = pairsAtoms2.size();

        if(ifSecondGasPresent){
            List<List<AtomType>> pairsAtomsGas2 = doPotential.getSpeciesPairs(speciesGas);
            int pairAtomSizeGas2 = pairsAtoms2.size();
        }

        box = this.makeBox();
        space = box.getSpace();
        box.getBoundary().setBoxSize(boxsize);
        if(ifMOPPresent){
            box.setNMolecules(speciesMOP,0);
        }

        List<Vector> oldPositions = new ArrayList<>();
        List<Vector> oldPositionsGas = new ArrayList<>();
        if(ifGraphenePresent) {
            box.addNewMolecule(speciesGrapheneOne);
            box.addNewMolecule(speciesGrapheneTwo);
            //box.addNewMolecule(speciesGrapheneThree);
            //box.addNewMolecule(speciesGrapheneFour);
        }
        box.setNMolecules(speciesGas, 1);
        if(ifSecondGasPresent)box.setNMolecules(speciesGasTwo, numGasTwo);

        IMolecule moleculeMOPOne = box.getMoleculeList().get(0);
        while (oldPositionsGas.size() < moleculeMOPOne.getChildList().size()) {
            oldPositionsGas.add(space.makeVector());
        }
        moleculeMOPOne.getChildList().forEach(atom ->{
            oldPositionsGas.get(atom.getIndex()).E(atom.getPosition());
            atom.getPosition().PE(centreMOP);
            Vector shift = box.getBoundary().centralImage(atom.getPosition());
            atom.getPosition().PE(shift);
        });


        //box.addNewMolecule(speciesGas);
      /*  Vector originSix = new Vector3D(0,-2, -2);
        IMolecule moleculeMOPOne = box.getMoleculeList().get(5);
        while (oldPositionsGas.size() < moleculeMOPOne.getChildList().size()) {
            oldPositionsGas.add(space.makeVector());
        }
        moleculeMOPOne.getChildList().forEach(atom ->{
            oldPositionsGas.get(atom.getIndex()).E(atom.getPosition());
            atom.getPosition().PE(originSix);
            Vector shift = box.getBoundary().centralImage(atom.getPosition());
            atom.getPosition().PE(shift);
        });
        Vector originSeven = new Vector3D(1,0, -1);
        IMolecule moleculeMOPTwo = box.getMoleculeList().get(6);
        moleculeMOPTwo.getChildList().forEach(atom ->{
            oldPositionsGas.get(atom.getIndex()).E(atom.getPosition());
            atom.getPosition().PE(originSeven);
            Vector shift = box.getBoundary().centralImage(atom.getPosition());
            atom.getPosition().PE(shift);
        });
        Vector originEight = new Vector3D(1,2, -2);
        IMolecule moleculeMOPThree= box.getMoleculeList().get(7);
        moleculeMOPThree.getChildList().forEach(atom ->{
            oldPositionsGas.get(atom.getIndex()).E(atom.getPosition());
            atom.getPosition().PE(originEight);
            Vector shift = box.getBoundary().centralImage(atom.getPosition());
            atom.getPosition().PE(shift);
        });*/
        box.getMoleculeList();
       // System.out.println(Constants.BOLTZMANN_K);
        //System.out.println(1/(Constants.BOLTZMANN_K*temperature));

        /*if(ifGraphenePresent && ifMOPPresent){
            if(ifSecondGasPresent){
                sm = new SpeciesManager.Builder().addSpecies(speciesMOP).addSpecies(speciesGrapheneOne).addSpecies(speciesGrapheneTwo).addSpecies(speciesGrapheneThree).addSpecies(speciesGrapheneFour).addSpecies(speciesGas).addSpecies(speciesGasTwo).build();
            }else {
                sm = new SpeciesManager.Builder().addSpecies(speciesMOP).addSpecies(speciesGrapheneOne).addSpecies(speciesGrapheneTwo).addSpecies(speciesGas).build();
            }
        } else if (ifGraphenePresent) {
            sm = new SpeciesManager.Builder().addSpecies(speciesGrapheneOne).addSpecies(speciesGrapheneTwo).addSpecies(speciesGas).build();
        } else {
            if(ifSecondGasPresent){
                sm = new SpeciesManager.Builder().addSpecies(speciesMOP).addSpecies(speciesGas).addSpecies(speciesGasTwo).build();
            }else {
                sm = new SpeciesManager.Builder().addSpecies(speciesMOP).addSpecies(speciesGas).build();
            }
        }*/
        if (ifGraphenePresent && ifMOPPresent){
            if (ifSecondGasPresent) {
                sm = new SpeciesManager.Builder().addSpecies(speciesMOP).addSpecies(speciesGrapheneOne).addSpecies(speciesGrapheneTwo).addSpecies(speciesGas).addSpecies(speciesGasTwo).build();
            }else {
                sm = new SpeciesManager.Builder().addSpecies(speciesMOP).addSpecies(speciesGrapheneOne).addSpecies(speciesGrapheneTwo).addSpecies(speciesGas).build();
            }
        } else if (ifGraphenePresent) {
            if (ifSecondGasPresent) {
                sm = new SpeciesManager.Builder().addSpecies(speciesGrapheneOne).addSpecies(speciesGrapheneTwo).addSpecies(speciesGas).addSpecies(speciesGasTwo).build();
            }else {
                sm = new SpeciesManager.Builder().addSpecies(speciesGrapheneOne).addSpecies(speciesGrapheneTwo).addSpecies(speciesGas).build();
            }
        } else if (ifMOPPresent) {
            if (ifSecondGasPresent) {
                sm = new SpeciesManager.Builder().addSpecies(speciesMOP).addSpecies(speciesGas).addSpecies(speciesGasTwo).build();
            }else {
                sm = new SpeciesManager.Builder().addSpecies(speciesMOP).addSpecies(speciesGas).build();
            }
        } else if (ifSecondGasPresent) {
            sm = new SpeciesManager.Builder().addSpecies(speciesGas).addSpecies(speciesGasTwo).build();
        }else {
            sm = new SpeciesManager.Builder().addSpecies(speciesGas).build();
        }
        UniversalSimulation.makeAtomPotentials(sm);
        PotentialMasterBonding pmBonding = new PotentialMasterBonding(sm, box);
        double nbrRange = truncatedRadiusLJ * 1.05 + 1;
        potentialMaster = new PotentialMasterList(getSpeciesManager(), box, 2, nbrRange, pmBonding.getBondingInfo());
        PotentialMasterCell potentialMasterCell = new PotentialMasterCell(getSpeciesManager(), box, 2, pmBonding.getBondingInfo());
        if(ifMOPPresent){
            ArrayList<ArrayList<Integer>> connectedAtoms1 = pdbReaderMOP.getConnectivity();
            //System.out.println(connectedAtoms1);
            ArrayList<ArrayList<Integer>> connectivityModified1 = pdbReaderMOP.getConnectivityModifiedWithoutRunning();
            //System.out.println(connectivityModified1);
            Map<Integer, String> atomMap1 = pdbReaderMOP.getAtomMapWithoutRunning();
            //System.out.println(atomMap1);
            // HashMap<Integer, String> atomMapModified1 = pdbReaderMOP.getAtomMapModifiedWithoutRunning();
            //System.out.println(atomMapModified1);
            ArrayList<Integer> bondList1 = pdbReaderMOP.getBondList(connectedAtoms1, atomMap1);
            // System.out.println(bondList1);
            Map<String, double[]> atomicPotMap1 = pdbReaderMOP.atomicPotMap();
            //System.out.println(atomicPotMap1);
            ArrayList<Integer> bondsNum1 = pdbReaderMOP.getBonds();

            //Map<Integer, String> atomIdentifierMapModified1 = PDBReader.atomIdentifierMapModified(connectivityModified1, atomMapModified1);
            Map<Integer, String> atomIdentifierMapModified1 = pdbReaderMOP.getModifiedAtomIdentifierMap();
            List<int[]> dupletsSorted1 = pdbReaderMOP.getDupletesSorted();
            List<int[]> tripletsSorted1 = pdbReaderMOP.getAnglesSorted();
            List<int[]> quadrupletsSorted1 = pdbReaderMOP.getTorsionSorted();

            Map<String[], List<int[]>> bondTypesMap1 = pdbReaderMOP.idenBondTypes(dupletsSorted1, atomIdentifierMapModified1);
            Map<String[], List<int[]>> angleTypesMap1 = pdbReaderMOP.idenAngleTypes(tripletsSorted1, atomIdentifierMapModified1);
            Map<String[], List<int[]>> torsionTypesMap1 = pdbReaderMOP.idenTorsionTypes(quadrupletsSorted1, atomIdentifierMapModified1);
            ArrayList<ArrayList<Integer>> modifiedOutput1 = new ArrayList<>();

            for (ArrayList<Integer> innerList : connectivityModified1) {
                ArrayList<Integer> modifiedInnerList = new ArrayList<>(innerList.subList(1, innerList.size()));
                modifiedOutput1.add(modifiedInnerList);
            }
            IntArrayList[] dupletsIntArrayList1 = new IntArrayList[modifiedOutput1.size()];

            for (i = 0; i < modifiedOutput1.size(); i++) {
                ArrayList<Integer> innerList = modifiedOutput1.get(i);
                IntArrayList intArrayList = new IntArrayList(innerList.size());
                for (int j = 0; j < innerList.size(); j++) {
                    intArrayList.add(innerList.get(j));
                }
                dupletsIntArrayList1[i] = intArrayList;
            }
            if(!ifMOPFixed){
                doPotential.setBondStretch(speciesMOP, bondTypesMap1, angleTypesMap1, torsionTypesMap1, bondsNum1, bondList1, quadrupletsSorted1, atomIdentifierMapModified1, atomicPotMap1, pmBonding);
            }
        }

        ArrayList<ArrayList<Integer>> connectedAtoms2 = pdbReaderReplicaNew.getConnectivity();
        ArrayList<ArrayList<Integer>> connectivityModified2 = pdbReaderReplicaNew.getConnectivityModifiedWithoutRunning();
        Map<Integer, String> atomMap2 = pdbReaderReplicaNew.getAtomMapWithoutRunning();
        ArrayList<Integer> bondList2 = pdbReaderReplicaNew.getBondList(connectedAtoms2, atomMap2);
        Map<String, double[]> atomicPotMap2 = pdbReaderReplicaNew.atomicPotMap();
        Map<Integer, String> atomIdentifierMapModified2 = pdbReaderReplicaNew.getatomIdentifierMapModified();
        List<int[]> dupletsSorted2 = pdbReaderReplicaNew.getDupletesSorted();
        List<int[]> tripletsSorted2 = pdbReaderReplicaNew.getAnglesSorted();
        List<int[]> quadrupletsSorted2 = pdbReaderReplicaNew.getTorsionSorted();
        ArrayList<Integer> bondsNum2 = pdbReaderReplicaNew.getBonds();

        Map<String[], List<int[]>> bondTypesMap2 = pdbReaderReplicaNew.idenBondTypes(dupletsSorted2, atomIdentifierMapModified2);
        Map<String[], List<int[]>> angleTypesMap2 = pdbReaderReplicaNew.idenAngleTypes(tripletsSorted2, atomIdentifierMapModified2);
        Map<String[], List<int[]>> torsionTypesMap2 = pdbReaderReplicaNew.idenTorsionTypes(quadrupletsSorted2, atomIdentifierMapModified2);

        ArrayList<ArrayList<Integer>> modifiedOutput2 = new ArrayList<>();
        for (ArrayList<Integer> innerList : connectivityModified2) {
            ArrayList<Integer> modifiedInnerList = new ArrayList<>(innerList.subList(1, innerList.size()));
            modifiedOutput2.add(modifiedInnerList);
        }
        IntArrayList[] dupletsIntArrayList2 = new IntArrayList[modifiedOutput2.size()];


        for (i = 0; i < modifiedOutput2.size(); i++) {
            ArrayList<Integer> innerList = modifiedOutput2.get(i);
            IntArrayList intArrayList = new IntArrayList(innerList.size());
            for (int j = 0; j < innerList.size(); j++) {
                intArrayList.add(innerList.get(j));
            }
            dupletsIntArrayList2[i] = intArrayList;
        }

        doPotential.setBondStretch(speciesGas, bondTypesMap2, angleTypesMap2, torsionTypesMap2, bondsNum2, bondList2, quadrupletsSorted2, atomIdentifierMapModified2, atomicPotMap2, pmBonding);

        if(ifSecondGasPresent){
            ArrayList<ArrayList<Integer>> connectedAtomsGas2 = pdbReaderReplicaGasNew.getConnectivity();
            ArrayList<ArrayList<Integer>> connectivityModifiedGas2 = pdbReaderReplicaGasNew.getConnectivityModifiedWithoutRunning();
            Map<Integer, String> atomMapGas2 = pdbReaderReplicaGasNew.getAtomMapWithoutRunning();
            ArrayList<Integer> bondListGas2 = pdbReaderReplicaGasNew.getBondList(connectedAtoms2, atomMap2);
            Map<String, double[]> atomicPotMapGas2 = pdbReaderReplicaGasNew.atomicPotMap();
            Map<Integer, String> atomIdentifierMapModifiedGas2 = pdbReaderReplicaGasNew.getatomIdentifierMapModified();
            List<int[]> dupletsSortedGas2 = pdbReaderReplicaGasNew.getDupletesSorted();
            List<int[]> tripletsSortedGas2 = pdbReaderReplicaGasNew.getAnglesSorted();
            List<int[]> quadrupletsSortedGas2 = pdbReaderReplicaGasNew.getTorsionSorted();
            ArrayList<Integer> bondsNumGas2 = pdbReaderReplicaGasNew.getBonds();

            Map<String[], List<int[]>> bondTypesMapGas2 = pdbReaderReplicaGasNew.idenBondTypes(dupletsSortedGas2, atomIdentifierMapModifiedGas2);
            Map<String[], List<int[]>> angleTypesMapGas2 = pdbReaderReplicaGasNew.idenAngleTypes(tripletsSortedGas2, atomIdentifierMapModifiedGas2);
            Map<String[], List<int[]>> torsionTypesMapGas2 = pdbReaderReplicaGasNew.idenTorsionTypes(quadrupletsSortedGas2, atomIdentifierMapModifiedGas2);

            ArrayList<ArrayList<Integer>> modifiedOutputGas2 = new ArrayList<>();
            for (ArrayList<Integer> innerList : connectivityModified2) {
                ArrayList<Integer> modifiedInnerList = new ArrayList<>(innerList.subList(1, innerList.size()));
                modifiedOutputGas2.add(modifiedInnerList);
            }
            IntArrayList[] dupletsIntArrayListGas2 = new IntArrayList[modifiedOutput2.size()];


            for (i = 0; i < modifiedOutputGas2.size(); i++) {
                ArrayList<Integer> innerList = modifiedOutputGas2.get(i);
                IntArrayList intArrayList = new IntArrayList(innerList.size());
                for (int j = 0; j < innerList.size(); j++) {
                    intArrayList.add(innerList.get(j));
                }
                dupletsIntArrayListGas2[i] = intArrayList;
            }

            if(ifSecondGasPresent)doPotential.setBondStretch(speciesGasTwo, bondTypesMapGas2, angleTypesMapGas2, torsionTypesMapGas2, bondsNumGas2, bondListGas2, quadrupletsSortedGas2, atomIdentifierMapModifiedGas2, atomicPotMapGas2, pmBonding);
            HashSet<AtomType> uniqueGasAtomTypes = new HashSet<AtomType>();
            if(ifSecondGasPresent )uniqueGasAtomTypes.addAll(speciesGasTwo.getAtomTypes());
            uniqueGasAtomTypes.addAll(speciesGas.getAtomTypes());
            System.out.println(uniqueGasAtomTypes);

            potentialMaster.doAllTruncationCorrection = false;
            System.out.println(speciesGas.getLeafAtomCount());
            List<List<AtomType>> pairsAtomsGas = doPotential.listFinal(uniqueGasAtomTypes);
            int uniqueGasAtomsSize = pairsAtomsGas.size();
            LJUFF[] p2LJGas = new LJUFF[uniqueGasAtomsSize];
            IPotential2[] p2ljGas = new IPotential2[uniqueGasAtomsSize];
            double[] sigmaIJGas = new double[uniqueGasAtomsSize];
            doPotential.doLJMD(pairsAtomsGas, p2LJGas, p2ljGas, truncatedRadiusLJ, potentialMaster, sigmaIJGas, potentialMasterCell);

            if(!ifGOMOPMove){
                HashSet<AtomType> uniqueAtomTypes = new HashSet<AtomType>();
                if(ifGraphenePresent) uniqueAtomTypes.addAll(speciesGrapheneOne.getAtomTypes());
                if(ifMOPPresent) uniqueAtomTypes.addAll(speciesMOP.getAtomTypes());
                System.out.println(uniqueAtomTypes);
                List<List<AtomType>> pairGasMembrane = doPotential.listFinal(uniqueGasAtomTypes, uniqueAtomTypes);
                int uniqueGasMembraneSize = pairGasMembrane.size();
                LJUFF[] p2LJGasMembrane = new LJUFF[uniqueGasMembraneSize];
                IPotential2[] p2ljGasMembrane = new IPotential2[uniqueGasMembraneSize];
                double[] sigmaIJGasMembrane = new double[uniqueGasMembraneSize];
                doPotential.doLJMD(pairGasMembrane, p2LJGasMembrane, p2ljGasMembrane, truncatedRadiusLJ, potentialMaster, sigmaIJGasMembrane, potentialMasterCell);
            }else {
                throw new RuntimeException("Error in setting up LJ as MOP or GO are moving");
            }

        }

/*

        if (ifGraphenePresent) {
            //GrapheneMOP
            //List<AtomType> listGrapheneIndividual = doPotential.listTwoSpeciesPairs(speciesGrapheneOne.getUniqueAtomTypes(), speciesGrapheneTwo.getUniqueAtomTypes());
            //pureGraphene
            //List<List<AtomType>> listGraphene = listMixedSpecial(speciesGrapheneOne.getUniqueAtomTypes(), speciesGrapheneTwo.getUniqueAtomTypes());
            //int pairAtomsGrapheneSize= listGraphene.size();
            // LJUFF[] p2LJGraphene = new LJUFF[pairAtomsGrapheneSize];
            // doLJ(listGraphene, potentialMasterCell, p2LJGraphene, pairAtomsGrapheneSize, truncatedRadiusLJ);

            //GrapheneGas
            List<List<AtomType>> listGrapheneGasFinal = doPotential.listGrapheneSpecial(speciesGrapheneOne.getUniqueAtomTypes(), speciesGas.getUniqueAtomTypes());
            //  List<List<AtomType>> listGrapheneGasFinal = SetPotential.listFinal(listGrapheneGasUnique);

            int pairAtomsGrapheneGasSize = listGrapheneGasFinal.size();
            LJUFF[] p2LJGrapheneGas = new LJUFF[pairAtomsGrapheneGasSize];
            IPotential2[] p2LJGraphenegas = new IPotential2[pairAtomsGrapheneGasSize];
            double[] sigmaGraphenegas = new double[pairAtomsGrapheneGasSize];
            P2Electrostatic[] P2Electrostatics = new P2Electrostatic[pairAtomsGrapheneGasSize];
            boolean doElectrostatics = false;
           // doPotential.doLJGrapheneMixed(listGrapheneGasFinal, p2LJGrapheneGas, p2LJGraphenegas,P2Electrostatics, truncatedRadiusLJ,potentialMaster, potentialMasterCell, sigmaGraphenegas, doElectrostatics);
        }

            List<AtomType> atomTypesMOP = null, atomTypeSecondGas = null, atomTypesgraphene = null,atomTypesGas;
            atomTypesGas = speciesGas.getUniqueAtomTypes();
            if(ifGraphenePresent) atomTypesgraphene = speciesGrapheneOne.getUniqueAtomTypes();
            if(ifMOPPresent) atomTypesMOP = speciesMOP.getAtomTypes();
            if(ifSecondGasPresent)atomTypeSecondGas = speciesGasTwo.getUniqueAtomTypes();
            // System.out.println(atomTypesGas + " " + atomTypesgraphene + " " +atomTypesMOP + " " +atomTypessecondGas);
            List<List<AtomType>> listSecondGas = new ArrayList<>();
            HashSet<AtomType> uniqueAtomTypes = new HashSet<AtomType>();
          /*  uniqueAtomTypes.addAll(atomTypesgraphene);
            uniqueAtomTypes.addAll(atomTypesGas);
            uniqueAtomTypes.addAll(atomTypesMOP);
            uniqueAtomTypes.addAll(atomTypeSecondGas);
            List<AtomType> uniqueList = new ArrayList<>(uniqueAtomTypes);

            if(ifGraphenePresent)doPotential.addUniqueElements(uniqueAtomTypes, atomTypesgraphene);
            if(ifMOPPresent)doPotential.addUniqueElements(uniqueAtomTypes, atomTypesMOP);
            doPotential.addUniqueElements(uniqueAtomTypes, atomTypesGas);
            if(ifSecondGasPresent) doPotential.addUniqueElements(uniqueAtomTypes, atomTypeSecondGas);

            // Convert HashSet back to List
            List<AtomType> uniqueList = new ArrayList<>(uniqueAtomTypes);
            // Display uniqueList
            System.out.println(uniqueList);
            listSecondGas = doPotential.listFinal(uniqueList);
            int pairAtomsGrapheneGasTwoSize = listSecondGas.size();
            LJUFF[] p2LJGrapheneGasTwo = new LJUFF[pairAtomsGrapheneGasTwoSize];
            IPotential2[] p2LJGraphenegasTwo = new IPotential2[pairAtomsGrapheneGasTwoSize];
            double[] sigmaGraphenegasTwo = new double[pairAtomsGrapheneGasTwoSize];
            P2Electrostatic[] P2ElectrostaticsTwo = new P2Electrostatic[pairAtomsGrapheneGasTwoSize];
            boolean doElectrostatics = false;
            doPotential.doLJGrapheneMixed(listSecondGas, p2LJGrapheneGasTwo, p2LJGraphenegasTwo,P2ElectrostaticsTwo, truncatedRadiusLJ,potentialMaster, potentialMasterCell, sigmaGraphenegasTwo, doElectrostatics);
*/
        box.getBoundary().setBoxSize(boxsize);
        final Vector bs = box.getBoundary().getBoxSize();
        System.out.println(bs);
        //integrator = new IntegratorVelocityVerlet(potentialMasterCell, random, 0.001, temperature, box);
        pcAgg = new PotentialComputeAggregate(pmBonding, potentialMaster);
        integratorMD = new IntegratorVelocityVerlet(pcAgg, random, 0.0005, temperature, box);
        integratorMD.setThermostatNoDrift(false);
        integratorMD.setIsothermal(false);
        nhc = new IntegratorListenerNHC(integratorMD, random, 3, 2);
        integratorMD.getEventManager().addListener(nhc);

        if (truncatedRadiusLJ > 0.5 * box.getBoundary().getBoxSize().getX(1)) {
            throw new RuntimeException("Truncation radius too large.  Max allowed is" + 0.5 * box.getBoundary().getBoxSize().getX(0));
        }
        potentialMasterCell.init();
      //  double u0 = potentialMasterCell.computeAll(false);

        double x = 1;
        Unit kcals = new UnitRatio(new PrefixedUnit(Prefix.KILO,Calorie.UNIT),Mole.UNIT);
     //   System.out.println( u0 + " "+ x +" inMain "  + kcals.fromSim(u0));
        pcAggMC = new PotentialComputeAggregate(pmBonding, potentialMaster);
        integratorMC = new IntegratorMC(pcAggMC, random, temperature, box);

        MCMoveMolecule translateMove = new MCMoveMolecule(random, pcAggMC, box);
        integratorMC.getMoveManager().addMCMove(translateMove);

        MCMoveMoleculeRotate rotateMove = new MCMoveMoleculeRotate(random, pcAggMC, box);
        integratorMC.getMoveManager().addMCMove(rotateMove);

         //System.exit(1);
      /*  while (u0 > 1e6*numMolecules) {
            System.out.println(x +" before");
            x *= 0.99;
             System.out.println( x +" =x");
            for (int j=0; j<pairAtomsMOPGasSize; j++){
                System.out.println(x + " " + sigmaIJ2[0]);
                p2LJMOPGas[j].setSigmaNew(x*sigmaIJ2[j]);
                ((P2SoftSphericalSumTruncatedForceShifted)p2LJMOPgas[j]).setTruncationRadius(truncatedRadiusLJ);
            }
            u0 = potentialMasterCell.computeAll(false);
            System.out.println( u0 + " inMain afterwards " +kcals.fromSim(u0));
        }
        integratorMC.reset();

        while (u0 > 1e4*numMolecules) {
            while (u0 > 1e4 * numMolecules) {
                integratorMC.doStep();
                u0 = integratorMC.getPotentialEnergy();
            }
            while (x < 1 && u0 <= 1e4 * numMolecules) {
                //  System.out.println("Inside Loop Two - 2");
                x /= 0.99;
                if (x > 1) x = 1;
                for(int j = 0; j< pairAtomsMOPGasSize; j++){
                    //  System.out.println("Inside Loop Two - 3");
                    p2LJMOPGas[j].setSigmaNew(x*sigmaMOPGas[j]);
                    ((P2SoftSphericalSumTruncatedForceShifted)p2LJMOPgas[j]).setTruncationRadius(truncatedRadiusLJ);
                }
                u0 = potentialMasterCell.computeAll(false);
                //System.out.println(u0 +" inside Array @");
            }

            integratorMC.reset();
        }*/
    }

    public static void main(String[] args) {
        MOPMDParams params = new MOPMDParams();
        ParseArgs.doParseArgs(params, args);
        int numAtomOne = params.numAtomOne;
        //int numAtomTwo = params.numAtomTwo;
        String confNameOne = params.confNameOne;
        String confNameGasOne = params.confNameGasOne;
        String confNameGasTwo = params.confNameGasTwo;
        String confNameGraphene = params.confNameGraphene;
        int numGasOne = params.numGasOne;
        int numGasTwo = params.numGasTwo;
        int truncatedRadiusCharge = params.truncatedRadius;
        int truncatedRadiusLJ = params.truncatedRadiusLJ;
        int temperature = params.temperature;
        Vector boxsize = params.boxSize;
        Vector grapheneOne = params.grapheneOne;
        Vector grapheneTwo = params.grapheneTwo;
        Vector grapheneThree = params.grapaheneThree;
        Vector grapheneFour = params.grapheneFour;
        Vector centreMOP = params.centreMOP;
        Vector centreMOPTwo = params.centreMOPTwo;
        boolean ifGraphenePresent = params.ifGraphenePresent;
        boolean ifGraphics = params.isGraphics;
        boolean isDynamic = params.isDynamic;
        boolean ifSecondGasPresent = params.ifSecondGasPresent;
        boolean setMOPmassInfinite = params.setMOPmassInfinite;
        double mu = params.mu;
        double sigma = params.sigma;
        int numMOP = params.numMOP;
        int numGas = params.numGas;
        int numGraphene = params.numGraphene;
        double density = params.density;
        int numSteps = params.numSteps;
        int numMolecules = numGas+numGraphene+numMOP;
        int mcSteps = params.mcSteps;
        boolean ifGrapheneFixed = params.ifGrapheneFixed;
        boolean ifMOPFixed = params.ifMOPFixed;
        boolean ifMOPPresent = params.ifMOPPresent;
        boolean ifGOMOPMove = params.ifGOMOPMove;
        MOPMD sim = new MOPMD(confNameOne, confNameGasOne,confNameGasTwo,numGasOne, numGasTwo, centreMOP, centreMOPTwo, confNameGraphene, grapheneOne, grapheneTwo,grapheneThree, grapheneFour,  temperature, truncatedRadiusCharge, truncatedRadiusLJ,sigma, mu, ifGraphenePresent, isDynamic, numMOP, numGraphene, numGas, boxsize, setMOPmassInfinite, ifGrapheneFixed, ifMOPPresent, ifMOPFixed, ifSecondGasPresent, ifGOMOPMove);
        MeterPotentialEnergyFromIntegrator meterU = new MeterPotentialEnergyFromIntegrator(sim.integratorMD);
        MeterTemperature meterT = new MeterTemperature(sim.box, 3);
        sim.potentialMaster.init();
        System.out.println("u0: "+sim.potentialMaster.computeAll(true)/numMolecules);
        MeterPressure meterP = new MeterPressure(sim.box, sim.pcAgg);
        meterP.setTemperature(temperature);
        double p0 = meterP.getDataAsScalar();
        System.out.println("p0: "+p0);

        DataProcessorForked dpZ = new DataProcessorForked() {
            DataDouble.DataInfoDouble dataInfo = new DataDouble.DataInfoDouble("Z", Null.DIMENSION);
            DataDouble data = new DataDouble();

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
        if (ifGraphics) {
            sim.getController().addActivity(new ActivityIntegrate(sim.integratorMD));
            final SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, "Octane MD", 3);
            DiameterHashByType dhbt = (DiameterHashByType) simGraphic.getDisplayBox(sim.box).getDiameterHash();
            if(ifMOPPresent){
                for (int i=0; i<sim.speciesMOP.getUniqueAtomTypes().size(); i++){
                    dhbt.setDiameter(sim.speciesMOP.getAtomType(i), 1);
                }
            }
           /* if(ifGraphenePresent){
                dhbt.setDiameter(sim.speciesGrapheneOne.getTypeByName("CX"), 1);
                dhbt.setDiameter(sim.speciesGrapheneOne.getTypeByName("CY"), 1);
                dhbt.setDiameter(sim.speciesGrapheneOne.getTypeByName("C4"), 1);
               // dhbt.setDiameter(sim.speciesGrapheneOne.getTypeByName("CZ"), 1);
                dhbt.setDiameter(sim.speciesGrapheneOne.getTypeByName("OJ"), 0.8);
                dhbt.setDiameter(sim.speciesGrapheneOne.getTypeByName("OL"), 0.8);
                dhbt.setDiameter(sim.speciesGrapheneOne.getTypeByName("OK"), 0.8);
                dhbt.setDiameter(sim.speciesGrapheneOne.getTypeByName("OE"), 0.8);
                dhbt.setDiameter(sim.speciesGrapheneOne.getTypeByName("HK"), 0.5);

              //  dhbt.setDiameter(sim.speciesGrapheneTwo.getAtomType(0), 0.5);
                ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGrapheneOne.getTypeByName("CX"), ColorExtra.darkGray);
                ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGrapheneOne.getTypeByName("CY"), ColorExtra.lightcyan);
                ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGrapheneOne.getTypeByName("C4"), ColorExtra.blue3);
               // ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGrapheneOne.getTypeByName("CZ"), ColorExtra.blue1);
                ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGrapheneOne.getTypeByName("OL"), ColorExtra.gold);
                ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGrapheneOne.getTypeByName("OJ"), ColorExtra.goldenRod);
                ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGrapheneOne.getTypeByName("OK"), ColorExtra.khaki);
                ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGrapheneOne.getTypeByName("OE"), ColorExtra.olive);
                ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGrapheneOne.getTypeByName("HK"), ColorExtra.tomato);
            }*/
          //  ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGas.getTypeByName("Ar"), ColorExtra.firebrick);
            ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGas.getAtomTypes().get(0), ColorExtra.SteelBlue);
            for (int i=0; i<sim.speciesGas.getUniqueAtomTypes().size(); i++){
                dhbt.setDiameter(sim.speciesGas.getAtomType(i), 1);
            }
            dhbt.setDiameter(sim.speciesGas.getAtomType(0), 2);
            if(ifGraphenePresent){
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
            }
            dhbt.setDiameter(sim.speciesMOP.getAtomType(0), 3.851);
            //dhbt.setDiameter(sim.speciesMOP.getAtomType(1),  3.5);
           // dhbt.setDiameter(sim.speciesMOP.getAtomType(2), 2.886);
           // dhbt.setDiameter(sim.speciesMOP.getAtomType(3), 3.851);
            //dhbt.setDiameter(sim.speciesMOP.getAtomType(4), 3.5);
            // dhbt.setDiameter(sim.speciesMOP.getAtomType(5), 2);
            // dhbt.setDiameter(sim.speciesMOP.getAtomType(5), 2);
            // dhbt.setDiameter(sim.speciesGas.getAtomType(0), 1.5);
            // dhbt.setDiameter(sim.speciesGas.getAtomType(1), 1);
            //  dhbt.setDiameter(sim.speciesSecondGas.getAtomType(0), 1.5);
            // dhbt.setDiameter(sim.speciesSecondGas.getAtomType(1), 1);

            //dhbt.setDiameter(sim.speciesGas.getAtomType(0), 1.2);
          //  dhbt.setDiameter(sim.speciesGas.getAtomType(1), 0.8);
           //  ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesMOP.getAtomType(1), ColorExtra.indianRed);
            ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesMOP.getAtomType(0), ColorExtra.copper);
            //((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesMOP.getAtomType(4), ColorExtra.red);
            //((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesMOP.getAtomType(2), ColorExtra.gray);
            //((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesMOP.getAtomType(3), ColorExtra.blue);
            // ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesMOP.getAtomType(2), ColorExtra.gray);

            ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGas.getAtomType(0), ColorExtra.lightGray);
           // ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGas.getAtomType(1), ColorExtra.lime);
            if(ifSecondGasPresent) ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGasTwo.getAtomType(0), ColorExtra.darkGray);
            //if(ifSecondGasPresent) ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGasTwo.getAtomType(1), ColorExtra.lightseagreen);
          // if(ifSecondGasPresent) ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGasTwo.getAtomType(1), ColorExtra.cornflowerblue);
            //((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGasTwo.getAtomType(1), ColorExtra.cornflowerblue);
            dhbt.setDiameter(sim.speciesGas.getAtomType(0), 3.85);
            //dhbt.setDiameter(sim.speciesGas.getAtomType(1), 3);
            if(ifSecondGasPresent)dhbt.setDiameter(sim.speciesGasTwo.getAtomType(0), 3.85);
            //dhbt.setDiameter(sim.speciesGasTwo.getAtomType(1),  2.886);

            simGraphic.getController().getReinitButton().setPostAction(simGraphic.getPaintAction(sim.box));

            DataSourceCountTime timer = new DataSourceCountTime(sim.integratorMD);
            DisplayTextBox timerBox = new DisplayTextBox();
            timerBox.setLabel("Time");
            DataPumpListener pumpSteps = new DataPumpListener(timer, timerBox, 100);
            sim.integratorMD.getEventManager().addListener(pumpSteps);
            simGraphic.add(timerBox);

            Unit perN = new SimpleUnit(Null.DIMENSION, numMolecules, "1/N", "1/N", false);

            AccumulatorHistory historyU = new AccumulatorHistory(new HistoryCollapsingDiscard());
            historyU.setTimeDataSource(timer);
            AccumulatorHistory historyU2 = new AccumulatorHistory(new HistoryCollapsingAverage());
            historyU2.setTimeDataSource(timer);
            AccumulatorAverageCollapsing avgEnergy = new AccumulatorAverageCollapsing();
            avgEnergy.setPushInterval(10);
            DataFork forkU = new DataFork(new IDataSink[]{historyU, historyU2, avgEnergy});
            DataPumpListener pumpU = new DataPumpListener(meterU, forkU, 10);
            sim.integratorMD.getEventManager().addListener(pumpU);
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

            if (sim.nhc != null) {
                IntegratorListenerNHC.DataSourceTotalEnergy meterTotalEnergy = new IntegratorListenerNHC.DataSourceTotalEnergy(sim.integratorMD, sim.nhc);
                AccumulatorHistory historyTotalEnergy = new AccumulatorHistory(new HistoryCollapsingDiscard());
                historyTotalEnergy.setTimeDataSource(timer);
                DataPumpListener pumpTotalEnergy = new DataPumpListener(meterTotalEnergy, historyTotalEnergy);
                sim.integratorMD.getEventManager().addListener(pumpTotalEnergy);
                DisplayPlotXChart plotTotalEnergy = new DisplayPlotXChart();
                plotTotalEnergy.setLabel("conserved energy");
                historyTotalEnergy.addDataSink(plotTotalEnergy.makeSink("total"));
                simGraphic.add(plotTotalEnergy);
            }

            AccumulatorHistory historyT = new AccumulatorHistory(new HistoryCollapsingDiscard());
            historyT.setTimeDataSource(timer);
            DataPumpListener pumpT = new DataPumpListener(meterT, historyT);
            sim.integratorMD.getEventManager().addListener(pumpT);
            DisplayPlotXChart plotT = new DisplayPlotXChart();
            plotT.setLabel("T");
            plotT.setUnit(Kelvin.UNIT);
            historyT.addDataSink(plotT.makeSink("T"));
            simGraphic.add(plotT);

            AccumulatorAverageCollapsing avgP = new AccumulatorAverageCollapsing();
            AccumulatorHistory historyP = new AccumulatorHistory(new HistoryCollapsingDiscard());
            historyP.setTimeDataSource(timer);
            AccumulatorHistory historyP2 = new AccumulatorHistory(new HistoryCollapsingAverage());
            historyP2.setTimeDataSource(timer);
            forkP.addDataSink(avgP);
            forkP.addDataSink(historyP);
            forkP.addDataSink(historyP2);
            DataPumpListener pumpP = new DataPumpListener(meterP, forkP, 10);
            sim.integratorMD.getEventManager().addListener(pumpP);
            DisplayPlotXChart plotP = new DisplayPlotXChart();
            plotP.setLabel("P");
            plotP.setUnit(Bar.UNIT);
            historyP.addDataSink(plotP.makeSink("P"));
            plotP.setLegend(new DataTag[]{historyP.getTag()}, "samples");
            historyP2.addDataSink(plotP.makeSink("Pavg"));
            plotP.setLegend(new DataTag[]{historyP2.getTag()}, "avg");
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

            simGraphic.makeAndDisplayFrame();
            return;
        }

        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integratorMD, numSteps/5));

        System.out.println("equilibration finished");
        long samples = numSteps / 100;
        long bs = samples / 100;
        if (bs == 0) bs = 1;

        AccumulatorAverageFixed accT = new AccumulatorAverageFixed(bs*4*10);
        DataPumpListener pumpT = new DataPumpListener(meterT, accT, 5);
        sim.integratorMD.getEventManager().addListener(pumpT);

        AccumulatorAverageFixed accU = new AccumulatorAverageFixed(bs*4*10);
        DataPumpListener pumpU = new DataPumpListener(meterU, accU, 5);
        sim.integratorMD.getEventManager().addListener(pumpU);

        AccumulatorAverageFixed accP = new AccumulatorAverageFixed(bs);
        forkP.addDataSink(accP);
        AccumulatorAverageFixed accZ = new AccumulatorAverageFixed(bs);
        dpZ.addDataSink(accZ);
        AccumulatorAverageFixed accZm1oR = new AccumulatorAverageFixed(bs);
        dpZm1oR.addDataSink(accZm1oR);
        DataPumpListener pumpP = new DataPumpListener(meterP, forkP, 20);
        sim.integratorMD.getEventManager().addListener(pumpP);

        int cycles = params.cycles;
        System.out.println(cycles+" cycles");
        long t1 = System.nanoTime();
        for (int i=0; i<cycles; i++) {
            sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integratorMD, numSteps/cycles));
            if (mcSteps>0) {
                sim.pcAggMC.init();
                sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integratorMD, mcSteps/cycles));
                sim.pcAgg.init();
            }
            double[] etaP = sim.nhc.getEtaP();
            sim.integratorMD.getEventManager().removeListener(sim.nhc);
            sim.nhc = new IntegratorListenerNHC(sim.integratorMD, sim.getRandom(), 3, 2);
            sim.nhc.setEtaP(etaP);
            sim.integratorMD.getEventManager().addListener(sim.nhc);
        }
        long t2 = System.nanoTime();

        IData dataT = accT.getData();
        double avgT = Kelvin.UNIT.fromSim(dataT.getValue(accT.AVERAGE.index));
        double errT = Kelvin.UNIT.fromSim(dataT.getValue(accT.ERROR.index));
        double corT = dataT.getValue(accT.BLOCK_CORRELATION.index);
        System.out.println("T: " + avgT + "   err: " + errT + "   cor: " + corT );
        double sdevT = Kelvin.UNIT.fromSim(dataT.getValue(accT.STANDARD_DEVIATION.index));
        System.out.println("sdev T: "+sdevT);

        IData dataU = accU.getData();
        double avgU = dataU.getValue(accU.AVERAGE.index) / numMolecules;
        double errU = dataU.getValue(accU.ERROR.index) / numMolecules;
        double corU = dataU.getValue(accU.BLOCK_CORRELATION.index);
        System.out.println("U: "+avgU+"   err: "+errU+"   cor: "+corU);

        IData dataP = accP.getData();
        double avgP = dataP.getValue(accP.AVERAGE.index);
        double errP = dataP.getValue(accP.ERROR.index);
        double corP = dataP.getValue(accP.BLOCK_CORRELATION.index);
        System.out.println("P: "+Bar.UNIT.fromSim(avgP)+"   err: "+errP+"   cor: "+corP);

        IData dataZ = accZ.getData();
        double avgZ = dataZ.getValue(accZ.AVERAGE.index);
        double errZ = dataZ.getValue(accZ.ERROR.index);
        double corZ = dataZ.getValue(accZ.BLOCK_CORRELATION.index);
        System.out.println("Z: "+avgZ+"   err: "+errZ+"   cor: "+corZ);

        IData dataZ_ = accZm1oR.getData();
        double avgZ_ = dataZ_.getValue(accZm1oR.AVERAGE.index);
        double errZ_ = dataZ_.getValue(accZm1oR.ERROR.index);
        double corZ_ = dataZ_.getValue(accZm1oR.BLOCK_CORRELATION.index);
        System.out.println("(Z-1)/rho: "+avgZ_+"   err: "+errZ_+"   cor: "+corZ_);

        System.out.println("time: "+(t2-t1)/1e9);
    }
    public Integrator getIntegrator(){
        return integratorMD;
    }

    public static class MOPMDParams extends ParameterBase {
        public double density = 0.003;
        public boolean setMOPmassInfinite = true;
        public int numGas = 1;
        public int numGraphene = 2;
        public int numMOP = 2;

        public boolean isDynamic = true;

        public int numAtomOne = 10;
        public int truncatedRadius = 15;
        public int temperature = 1200;
        public double sigma = 5.2;
        public int numAtomTwo = 1;
        public int numSteps = 10000000;
        public double mu = -2500;
        public int mcSteps=0;
        public int cycles=5000;
        public int truncatedRadiusLJ = 15;
        public String confNameOne = "F://Avagadro//mop//tetra_cu";
        public String confNameGasTwo = "F://Avagadro//molecule//He";
        public String confNameGraphene = "F://Avagadro//holeGO60";
      // public String confNameGraphene = "F://Avagadro//molecule//ethaneholeGO60";
        public String confNameGasOne = "F://Avagadro//molecule//Ar";
        public Vector centreMOPTwo = new Vector3D(50,50,50);
        public Vector centreMOP = new Vector3D(-5,-5,5);
        public Vector boxSize = new Vector3D(30,30,50);
        public Vector grapheneOne = new Vector3D(15.0,15.0,10.0);
        public Vector grapheneTwo = new Vector3D(15.0,15.0, -10);
        public Vector grapaheneThree = new Vector3D(15.0,15.0,25);
        public Vector grapheneFour = new Vector3D(15.0,15.0,-13);
        public boolean isGraphics = true;
        public boolean ifGraphenePresent = true;
        public boolean ifGrapheneFixed = true;
        public boolean ifSecondGasPresent = true;
        public boolean ifMOPPresent = true;
        public boolean ifGOMOPMove = false;
        public boolean ifMOPFixed = true;
        public int numGasOne =2 ;
        public int numGasTwo =0;
    }
}