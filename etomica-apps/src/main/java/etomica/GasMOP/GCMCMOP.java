package etomica.GasMOP;

import etomica.action.BoxInflate;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.DiameterHashByType;
import etomica.box.Box;
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
import etomica.molecule.MoleculeSourceRandomMolecule;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.potential.*;
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

import java.awt.*;
import java.util.*;
import java.util.List;

public class GCMCMOP extends Simulation {
    public IntegratorMC integrator;
    public MCMoveMolecule mcMoveMolecule;
    public MCMoveMoleculeRotate mcMoveMoleculeRotate;
    public MCMoveInsertDelete mcMoveID, mcMoveIDTwo;
    public ISpecies speciesMOP, speciesGas, speciesGrapheneOne, speciesGrapheneTwo, speciesGrapheneThree, speciesGrapheneFour,speciesGrapheneFive,speciesGrapheneSix,speciesGrapheneSeven,speciesGrapheneEight, speciesSecondGas,speciesGrapheneNine,speciesGrapheneTen,speciesGrapheneEleven,speciesGrapheneTwelve, speciesGrapheneThirteen;
    public Box box;
    public P2LennardJones potential;
    public SpeciesManager sm;
    public static int atom1, atom2, atom3;
    public static double sigmaIKey, sigmaJKey, epsilonIKey, epsilonJKey;
    public int i;
    public static String atomName1,atomName2, atomName3 ;
    public static List<List<AtomType>> listMOPGasMixed = new ArrayList<>();
    public static List<List<AtomType>> listGrapheneGasMixed = new ArrayList<>();
    public GCMCMOP(String confNameOne, String confNameGasOne, String confNameGasTwo, String confNameGraphene, Vector grapheneOne, Vector grapheneTwo,Vector grapheneThree,Vector grapheneFour,Vector grapheneFive,Vector grapheneSix,Vector grapheneSeven,Vector grapheneEight,Vector grapheneNine, Vector grapheneTen, Vector grapheneEleven, Vector grapheneTwelve, Vector grapheneThirteen, int numMolOne, int temperature, int truncatedRadius, int truncatedRadiusLJ, double sigma, double mu1, double mu2, boolean ifGraphenePresent, boolean ifSecondGasPresent, boolean ifMultipleGraphenePresent, boolean ifMoveRotateMoves, Vector boxSize){
        super(Space3D.getInstance());

        //Make Species
        setRandom(new RandomMersenneTwister(11));
        PDBReaderMOP pdbReaderMOP = new PDBReaderMOP();
        PDBReaderReplica pdbReaderReplica = new PDBReaderReplica();
        PDBReaderReplicaNew pdbReaderReplicaNew = new PDBReaderReplicaNew();
        grapheneReader grapheneReader = new grapheneReader();
        DoPotential doPotential = new DoPotential();
        /*speciesMOP = pdbReaderMOP.getSpecies(confNameOne, false);
        addSpecies(speciesMOP);*/
        if(ifGraphenePresent){
            speciesGrapheneOne = grapheneReader.getSpecies(confNameGraphene, grapheneOne);
            speciesGrapheneFive = grapheneReader.getSpecies(confNameGraphene, grapheneFive);
            addSpecies(speciesGrapheneOne);
            addSpecies(speciesGrapheneFive);
            if(ifMultipleGraphenePresent){
                speciesGrapheneTwo = grapheneReader.getSpecies(confNameGraphene, grapheneTwo);
                speciesGrapheneThree = grapheneReader.getSpecies(confNameGraphene, grapheneThree);
                speciesGrapheneFour = grapheneReader.getSpecies(confNameGraphene, grapheneFour);
                speciesGrapheneSix = grapheneReader.getSpecies(confNameGraphene, grapheneSix);
                speciesGrapheneSeven = grapheneReader.getSpecies(confNameGraphene, grapheneSeven);
                speciesGrapheneEight = grapheneReader.getSpecies(confNameGraphene, grapheneEight);
                speciesGrapheneNine = grapheneReader.getSpecies(confNameGraphene, grapheneNine);
                speciesGrapheneTen = grapheneReader.getSpecies(confNameGraphene, grapheneTen);
                speciesGrapheneEleven = grapheneReader.getSpecies(confNameGraphene, grapheneEleven);
                speciesGrapheneTwelve = grapheneReader.getSpecies(confNameGraphene, grapheneTwelve);
                speciesGrapheneThirteen = grapheneReader.getSpecies(confNameGraphene, grapheneThirteen);
                addSpecies(speciesGrapheneTwo);
                addSpecies(speciesGrapheneThree);
                addSpecies(speciesGrapheneFour);
                addSpecies(speciesGrapheneSix);
                addSpecies(speciesGrapheneSeven);
                addSpecies(speciesGrapheneEight);
                addSpecies(speciesGrapheneNine);
                addSpecies(speciesGrapheneTen);
                addSpecies(speciesGrapheneEleven);
                addSpecies(speciesGrapheneTwelve);
                addSpecies(speciesGrapheneThirteen);
            }
        }
        //species Gas
        speciesGas = pdbReaderReplica.getSpecies(confNameGasOne);
        addSpecies(speciesGas);

        //SecondGas and Bonding Parameters of Second Gas
       if(ifSecondGasPresent){
            speciesSecondGas = pdbReaderReplicaNew.getSpecies(confNameGasTwo);
            addSpecies(speciesSecondGas);
            pairsAtomsSecondGas = doPotential.getSpeciesPairs(speciesSecondGas);
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

       //Make Box and add molecules
        box = this.makeBox();
       // box.setNMolecules(speciesMOP, numMolOne);
        if(ifGraphenePresent ){
            if(ifMultipleGraphenePresent ){
                box.addNewMolecule(speciesGrapheneOne);
                box.addNewMolecule(speciesGrapheneFive);
                box.addNewMolecule(speciesGrapheneTwo);
                box.addNewMolecule(speciesGrapheneThree);
                box.addNewMolecule(speciesGrapheneFour);
                box.addNewMolecule(speciesGrapheneSix);
                box.addNewMolecule(speciesGrapheneSeven);
                box.addNewMolecule(speciesGrapheneEight);
                box.addNewMolecule(speciesGrapheneNine);
                box.addNewMolecule(speciesGrapheneTen);
                box.addNewMolecule(speciesGrapheneEleven);
                box.addNewMolecule(speciesGrapheneTwelve);
                box.addNewMolecule(speciesGrapheneThirteen);
            }else {
                box.addNewMolecule(speciesGrapheneOne);
                box.addNewMolecule(speciesGrapheneFive);
            }
        }
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
        //sm = new SpeciesManager.Builder().addSpecies(speciesGrapheneOne).addSpecies(speciesGrapheneFive).addSpecies(speciesGrapheneTwo).addSpecies(speciesGrapheneThree).addSpecies(speciesGrapheneFour). addSpecies(speciesGrapheneSix). addSpecies(speciesGrapheneSeven). addSpecies(speciesGrapheneEight).addSpecies(speciesGrapheneThirteen).addSpecies(speciesGas).addSpecies(speciesSecondGas).build();
        sm = new SpeciesManager.Builder().addSpecies(speciesGrapheneOne).addSpecies(speciesGrapheneFive).addSpecies(speciesGrapheneTwo).addSpecies(speciesGrapheneThree).addSpecies(speciesGrapheneFour). addSpecies(speciesGrapheneSix). addSpecies(speciesGrapheneSeven). addSpecies(speciesGrapheneEight).addSpecies(speciesGrapheneNine).addSpecies(speciesGrapheneTen).addSpecies(speciesGrapheneEleven).addSpecies(speciesGrapheneTwelve).addSpecies(speciesGrapheneThirteen).addSpecies(speciesGas).addSpecies(speciesSecondGas).build();
        //Bonding parameters MOP
        /*ArrayList<ArrayList<Integer>> connectedAtoms1 = pdbReaderMOP.getConnectivityWithoutRunning();
        ArrayList<ArrayList<Integer>> connectivityModified1 = pdbReaderMOP.getConnectivityModifiedWithoutRunning();
        Map<Integer,String> atomMap1 = pdbReaderMOP.getAtomMapWithoutRunning();
        HashMap<Integer, String> atomMapModified1 = pdbReaderMOP.getAtomMapModifiedWithoutRunning();
        ArrayList<Integer> bondList1 = pdbReaderMOP.getBondList(connectedAtoms1, atomMap1);
        Map<String, double[]> atomicPotMap1 = pdbReaderMOP.atomicPotMap();
        ArrayList<Integer> bondsNum1 = pdbReaderMOP.getBonds();
        Map<Integer, String> atomIdentifierMapModified1 = pdbReaderMOP.getModifiedAtomIdentifierMap();
        List<int[]> dupletsSorted1= pdbReaderMOP.getDupletesSorted();
        List<int[]>tripletsSorted1= pdbReaderMOP.getAnglesSorted();
        List<int[]>quadrupletsSorted1= pdbReaderMOP.getTorsionSorted();
        Map<String[],List<int[]>> bondTypesMap1= pdbReaderMOP.idenBondTypes(dupletsSorted1, atomIdentifierMapModified1);
        Map<String[],List<int[]>> angleTypesMap1= pdbReaderMOP.idenAngleTypes(tripletsSorted1, atomIdentifierMapModified1);
        Map<String[],List<int[]>> torsionTypesMap1= pdbReaderMOP.idenTorsionTypes(quadrupletsSorted1, atomIdentifierMapModified1);*/

        //Bonding Parameters GasOne
        ArrayList<ArrayList<Integer>> connectedAtoms2 =pdbReaderReplica.getConnectivityWithoutRunning();
        ArrayList<ArrayList<Integer>> connectivityModified2 = pdbReaderReplica.getConnectivityModifiedWithoutRunning();
        Map<Integer,String> atomMap2 = pdbReaderReplica.getAtomMapWithoutRunning();
        ArrayList<Integer> bondList2 = pdbReaderReplica.getBondList(connectedAtoms2, atomMap2);
        Map<String, double[]> atomicPotMap2 = pdbReaderReplica.atomicPotMap();
        Map<Integer, String> atomIdentifierMapModified2 = pdbReaderReplica.getatomIdentifierMapModified();
        List<int[]>dupletsSorted2= pdbReaderReplica.getDupletesSorted();
        List<int[]>tripletsSorted2=pdbReaderReplica.getAnglesSorted();
        List<int[]>quadrupletsSorted2=pdbReaderReplica.getTorsionSorted();
        ArrayList<Integer> bondsNum2 = pdbReaderReplica.getBonds();
        Map<String[],List<int[]>> bondTypesMap2= pdbReaderReplica.idenBondTypes(dupletsSorted2, atomIdentifierMapModified2);
        Map<String[],List<int[]>> angleTypesMap2= pdbReaderReplica.idenAngleTypes(tripletsSorted2, atomIdentifierMapModified2);
        Map<String[],List<int[]>> torsionTypesMap2= pdbReaderReplica.idenTorsionTypes(quadrupletsSorted2, atomIdentifierMapModified2);

        //SetBonding Potential
        UniversalSimulation.makeAtomPotentials(sm);
        PotentialMasterBonding pmBonding = new PotentialMasterBonding(sm, box);
        PotentialMasterCell potentialMasterCell = new PotentialMasterCell(getSpeciesManager(), box, 2, pmBonding.getBondingInfo());
        //MOP
       // doPotential.doBondStrech(speciesMOP, bondTypesMap1, angleTypesMap1, torsionTypesMap1,bondsNum1,bondList1, quadrupletsSorted1, atomIdentifierMapModified1,atomicPotMap1, pmBonding);
        //GasOne
        doPotential.doBondStrech(speciesGas,bondTypesMap2, angleTypesMap2, torsionTypesMap2, bondsNum2, bondList2,quadrupletsSorted2, atomIdentifierMapModified2, atomicPotMap2, pmBonding);
        //GasTwo
        if(ifSecondGasPresent){
            doPotential.doBondStrech(speciesSecondGas,bondTypesMapSecondGas, angleTypesMapSecondGas, torsionTypesMapSecondGas, bondsNumSecondGas, bondListSecondGas,quadrupletsSortedSecondGas, atomIdentifierMapModifiedSecondGas, atomicPotMapSecondGas, pmBonding);
        }
        List<AtomType> listGas;

        //NonBonded
        potentialMasterCell.doAllTruncationCorrection = false;
        if(ifSecondGasPresent){
            listGas = doPotential.listTwoSpeciesPairs(speciesGas.getUniqueAtomTypes(), speciesSecondGas.getUniqueAtomTypes());
        } else {
            listGas = speciesGas.getUniqueAtomTypes();
        }
        List<AtomType> listMOPGas = doPotential.listTwoSpeciesPairs(listGas, listGas);
        //MOP-Gas
      //  List<AtomType> listMOPGas = doPotential.listTwoSpeciesPairs(speciesMOP.getUniqueAtomTypes(), listGas);
        System.out.println(listMOPGas);
        List<List<AtomType>> listMOPGasPairs = doPotential.listFinal(listMOPGas);
        LJUFF[] p2LJMOPGas = new LJUFF[listMOPGasPairs.size()];
        doPotential.doLJ(listMOPGasPairs, potentialMasterCell,  p2LJMOPGas, listMOPGasPairs.size(), truncatedRadiusLJ);

        //Graphene-Gas
        List<List<AtomType>> listGrapheneGas = doPotential.listGrapheneSpecial(speciesGrapheneOne.getUniqueAtomTypes(), listGas);
        LJUFF[] p2LJGrapheneGas = new LJUFF[listGrapheneGas.size()];
        doPotential.doLJ(listGrapheneGas, potentialMasterCell,  p2LJGrapheneGas, listGrapheneGas.size(), truncatedRadiusLJ);

        integrator = new IntegratorMC(potentialMasterCell, random, temperature, box);
        potentialMasterCell.doOneTruncationCorrection = true;
        if(ifMoveRotateMoves){
            mcMoveMolecule = new MCMoveMolecule(random, potentialMasterCell, box);
            ((MoleculeSourceRandomMolecule) mcMoveMolecule.getMoleculeSource()).setSpecies(speciesGas);
            mcMoveMolecule.setStepSize(0.2 * sigma);
            ((MCMoveStepTracker) mcMoveMolecule.getTracker()).setTunable(false);
            integrator.getMoveManager().addMCMove(mcMoveMolecule);
            mcMoveMoleculeRotate = new MCMoveMoleculeRotate(random, potentialMasterCell, box);
            ((MoleculeSourceRandomMolecule) mcMoveMoleculeRotate.getMoleculeSource()).setSpecies(speciesGas);
            mcMoveMoleculeRotate.setStepSize(0.2 * sigma);
            ((MCMoveStepTracker) mcMoveMoleculeRotate.getTracker()).setTunable(false);
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

        box.getBoundary().setBoxSize(boxSize);
        final Vector bs = box.getBoundary().getBoxSize();
        System.out.println(bs);
        integrator.getMoveManager().setEquilibrating(false);

        BoxInflate inflater = new BoxInflate(box, space);
        inflater.actionPerformed();
        potential = new P2LennardJones(sigma, 1.0);
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
        int truncatedRadius = params.truncatedRadius;
        int truncatedRadiusLJ = params.truncatedRadiusLJ;
        int temperature = params.temperature;
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
        boolean ifGraphenePresent = params.ifGraphenePresent;
        boolean ifSecondGasPresent = params.ifSecondGasPresent;
        boolean ifMultipleGraphenePresent = params.ifMultipleGraphenePresent;
        boolean ifMoveRotateMoves = params.ifMoveRotateMoves;
        boolean doGraphics = params.doGraphics;
        double mu1 = params.mu1;
        double mu2= params.mu2;
        double sigma = params.sigma;
        GCMCMOP sim = new GCMCMOP(confNameOne, confNameGasOne, confNameGasTwo, confNameGraphene, grapheneOne, grapheneTwo, grapheneThree, grapheneFour, grapheneFive, grapheneSix, grapheneSeven, grapheneEight,grapheneNine, grapheneTen, grapheneEleven, grapheneTwelve, grapheneThirteen,numAtomOne,  temperature, truncatedRadius, truncatedRadiusLJ,sigma, mu1,mu2, ifGraphenePresent, ifSecondGasPresent,ifMultipleGraphenePresent,ifMoveRotateMoves, boxSize);
        if(doGraphics){
            final SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, "GCMC");
            DiameterHashByType dhbt = (DiameterHashByType) simGraphic.getDisplayBox(sim.box).getDiameterHash();
           /* dhbt.setDiameter(sim.speciesMOP.getAtomType(0), 1);
            dhbt.setDiameter(sim.speciesMOP.getAtomType(1), 1);
            dhbt.setDiameter(sim.speciesMOP.getAtomType(2), 1);
            dhbt.setDiameter(sim.speciesMOP.getAtomType(3), 1);
            dhbt.setDiameter(sim.speciesMOP.getAtomType(4), 1);
            dhbt.setDiameter(sim.speciesMOP.getAtomType(5), 1);*/
            if(ifGraphenePresent){
                dhbt.setDiameter(sim.speciesGrapheneOne.getAtomType(0), 1.5);
                /*dhbt.setDiameter(sim.speciesGrapheneOne.getAtomType(1), 1.5);
                dhbt.setDiameter(sim.speciesGrapheneOne.getAtomType(2), 1.5);
                dhbt.setDiameter(sim.speciesGrapheneOne.getAtomType(3), 1);
                dhbt.setDiameter(sim.speciesGrapheneOne.getAtomType(4), 1.5);
                dhbt.setDiameter(sim.speciesGrapheneOne.getAtomType(5), 1.5);
                dhbt.setDiameter(sim.speciesGrapheneOne.getAtomType(6), 1.5);
                dhbt.setDiameter(sim.speciesGrapheneOne.getAtomType(7), 1.5);*/
                dhbt.setDiameter(sim.speciesGas.getAtomType(0),2);
                dhbt.setDiameter(sim.speciesGas.getAtomType(0),1);
                ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGrapheneOne.getTypeByName("CX"), Color.lightGray);
                ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGrapheneFive.getTypeByName("CX"), Color.lightGray);
           }
            if (ifGraphenePresent){
                ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGrapheneOne.getAtomType(0), Color.gray);
              /* ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGrapheneOne.getAtomType(1), Color.gray);
                ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGrapheneOne.getAtomType(2), ColorExtra.lightcyan);
                ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGrapheneOne.getAtomType(3), Color.red);
                ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGrapheneOne.getAtomType(4), Color.gray);
                ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGrapheneOne.getAtomType(5), ColorExtra.lightcyan);
                ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGrapheneOne.getAtomType(6), ColorExtra.lightcyan);
                ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGrapheneOne.getAtomType(7), ColorExtra.lightcyan);*/
            }
            dhbt.setDiameter(sim.speciesGas.getAtomType(0),2);
            dhbt.setDiameter(sim.speciesGas.getAtomType(1),1);
           // dhbt.setDiameter(sim.speciesSecondGas.getAtomType(0),2);

            //MOP
         /*   dhbt.setDiameter(sim.speciesMOP.getAtomType(0), 5);
            dhbt.setDiameter(sim.speciesMOP.getAtomType(1), 2);
            dhbt.setDiameter(sim.speciesMOP.getAtomType(2), 2);
            dhbt.setDiameter(sim.speciesMOP.getAtomType(3), 2);
            dhbt.setDiameter(sim.speciesMOP.getAtomType(4), 2);
            dhbt.setDiameter(sim.speciesMOP.getAtomType(4), 2);
            ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesMOP.getAtomType(0), ColorExtra.goldenRod);
            ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesMOP.getAtomType(1), ColorExtra.lightcyan);
            ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesMOP.getAtomType(2), ColorExtra.darkRed);
            ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesMOP.getAtomType(3), ColorExtra.gray);
            ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesMOP.getAtomType(4), ColorExtra.lightcyan);
            ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesMOP.getAtomType(5), ColorExtra.gray);*/

           //((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesMOP.getTypeByName("Cu"), Color.white);
            //((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesMOP.getTypeByName("H"), Color.RED);
          // ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGas.getTypeByName("C_2"), Color.cyan);
           // ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGas.getTypeByName("O_2"), Color.blue);

           ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGas.getAtomType(0), ColorExtra.cornflowerblue);
           ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesGas.getAtomType(1), Color.cyan);
            //((ColorSchemeByType) simGraphic.getDisplayBox(sim.box()).getColorScheme()).setColor(sim.speciesSecondGas.getAtomType(0), ColorExtra.seaGreen);

            simGraphic.makeAndDisplayFrame("GCMC");
            ActivityIntegrate ai2 = new ActivityIntegrate(sim.integrator);
            sim.getController().addActivity(ai2, Long.MAX_VALUE, 1.0);
            return;
        }
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, params.numSteps / 10));

        int pInterval = 2 ;
        int bs = params.numSteps / (pInterval * 50);
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
        DataPumpListener pumpDensity = new DataPumpListener(densityMeter, densityAccumulator);
        sim.integrator.getEventManager().addListener(pumpDensity);

        long t1 = System.currentTimeMillis();
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, params.numSteps));
        long t2 = System.currentTimeMillis();
        System.out.println("runtime: " + (t2 - t1) * 0.001);
        Unit pUnit = Bar.UNIT;
        double avgP = pAccumulator.getData(pAccumulator.AVERAGE).getValue(0);
        double errP = pAccumulator.getData(pAccumulator.ERROR).getValue(0);
        double corP = pAccumulator.getData(pAccumulator.BLOCK_CORRELATION).getValue(0);
        System.out.println("P " + pUnit.fromSim(avgP) + " " + pUnit.fromSim(errP) + " " + pUnit.fromSim(corP));

        double avgRho = densityAccumulator.getData(densityAccumulator.AVERAGE).getValue(0);
        double errRho = densityAccumulator.getData(densityAccumulator.ERROR).getValue(0);
        double corRho = densityAccumulator.getData(densityAccumulator.BLOCK_CORRELATION).getValue(0);
        Unit m3 = new CompoundUnit(new Unit[] {Meter.UNIT}, new double[] {3});
        Unit kgm3 = new UnitRatio(new PrefixedUnit(Prefix.KILO, Gram.UNIT), m3);
        double numAtomsAvg = avgRho * sim.box.getBoundary().volume();
        System.out.println( "Num atoms : " +numAtomsAvg);

        double avgPE = energyAccumulator.getData(energyAccumulator.AVERAGE).getValue(0) / numAtomsAvg;
        double errPE = energyAccumulator.getData(energyAccumulator.ERROR).getValue(0) / numAtomsAvg;
        double corPE = energyAccumulator.getData(energyAccumulator.BLOCK_CORRELATION).getValue(0);
        System.out.println("PE " + avgPE + " " + errPE + " " + corPE);


        double numAtomGasOne = nMoleculesGasOne.getDataAsScalar();
        System.out.println("Atoms gas One: " +numAtomGasOne);

        if(ifSecondGasPresent){
            double numAtomGasTwo = nMoleculesGasTwo.getDataAsScalar();
            System.out.println("Atoms gas Two: " +numAtomGasTwo);
        }



        System.out.println("rho " + kgm3.fromSim(avgRho) + " " + errRho + " " + corRho);

        // expected values based on 10^8 steps
        // stdev based on 10^8 steps for N=500 uncertainty, scaled up to 10^6 steps
        //   other sims have significant correlation
        // 4 sigma should fail 1 in 16,000 runs

        double expectedP = 0.0815; // finite size effect smaller than uncertainty
        double stdevP = 0.03;
        if (Double.isNaN(avgP) || Math.abs(avgP - expectedP) / stdevP > 4) {
            System.exit(1);
        }

        double expectedPE = -4.492; // finite size effect comparable to uncertainty
        double stdevPE = 0.09;
        if (Double.isNaN(avgPE) || Math.abs(avgPE - expectedPE) / stdevPE > 4) {
            System.exit(2);
        }

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
        public String confNameOne = "F://Avagadro//mop//tetra_cu";
        public String confNameGasOne = "F://Avagadro//molecule//ch4";
        public String confNameGasTwo ="F://Avagadro//molecule//ethane";
        public String confNameGraphene = "F://Avagadro//pristine";
        public Vector grapheneThirteen = new Vector3D(15.0,15.0, 0.0);
        public Vector grapheneOne = new Vector3D(15.0,15.0, 10.0);
        public Vector grapheneTwo = new Vector3D(15.0,15.0, 20.0);
        public Vector grapheneThree = new Vector3D(15.0,15.0, 30.0);
        public Vector grapheneFour = new Vector3D(15.0,15.0, 40.0);
        public Vector grapheneNine = new Vector3D(15.0,15.0, 50.0);
        public Vector grapheneTen = new Vector3D(15.0,15.0, 60.0);
        public Vector grapheneFive = new Vector3D(15.0,15.0, -10.0);
        public Vector grapheneSix = new Vector3D(15.0,15.0, -20.0);
        public Vector grapheneSeven = new Vector3D(15.0,15.0, -30.0);
        public Vector grapheneEight = new Vector3D(15.0,15.0, -40.0);
        public Vector grapheneEleven = new Vector3D(15.0,15.0, -50.0);
        public Vector grapheneTwelve = new Vector3D(15.0,15.0, -60.0);

        public Vector boxSize = new Vector3D(35.0,32.0,140);
        public boolean ifGraphenePresent = true;
        public boolean ifMultipleGraphenePresent = true;
        public boolean ifSecondGasPresent = true;
        public boolean ifMoveRotateMoves = true;
        public boolean doGraphics = true;
        public int numAtomOne = 1;
        public int truncatedRadius = 5;
        public int temperature = 270;
        public double sigma = 5.2;
        public int numAtomTwo = 1;
        public int numSteps = 10000;
        public double mu1 = -2500;
        public double mu2 = -3000;
        public int truncatedRadiusLJ = 10;
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
