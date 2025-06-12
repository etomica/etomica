/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;

import etomica.action.PDBWriter;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.chem.elements.*;
import etomica.graph.model.Graph;
import etomica.graphics.ColorSchemeRandomByMolecule;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayBoxCanvasG3DSys;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorEvent;
import etomica.integrator.IntegratorListener;
import etomica.integrator.mcmove.MCMoveStepTracker;
import etomica.math.SpecialFunctions;
import etomica.molecule.IMoleculeList;
import etomica.molecule.MoleculePositionCOM;
import etomica.potential.*;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesBuilder;
import etomica.species.SpeciesManager;
import etomica.units.Degree;
import etomica.units.Electron;
import etomica.units.Kelvin;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.util.collections.IntArrayList;
import etomica.util.random.RandomMersenneTwister;
import etomica.virial.MayerFunction;
import etomica.virial.MayerHardSphere;
import etomica.virial.MayerMix;
import etomica.virial.cluster.*;
import etomica.virial.mcmove.*;
import etomica.virial.wheatley.ClusterWheatleyHS;
import etomica.virial.wheatley.ClusterWheatleySoftDerivatives;
import etomica.virial.wheatley.ClusterWheatleySoftDerivativesBD;

import java.awt.*;
import java.util.Arrays;
import java.util.Collections;
import java.util.Map;
import java.util.Set;

/**
 * Compute pure, binary, ternary and quaternary mixture virial coefficients using overlap sampling simulations
 * for some molecules using the TraPPE force fields.
 *
 */
public class VirialTraPPE {


    public static void main(String[] args) {
        VirialParam params = new VirialParam();
        boolean isCommandline = args.length > 0;
        if (isCommandline) {
            ParseArgs.doParseArgs(params, args);
        } else {
            // Customize Interactive Parameters Here
            params.chemForm = new ChemForm[]{ChemForm.isobutanol};
            params.nPoints = 2; //B order
            params.temperature = 450;
            params.diagram = "BC";
            params.numSteps = 2000000000;
            params.refFrac = -1;
            params.seed = new int[]{-997863083, -271031748, 477857854, 450871245};
            params.dorefpref = false;
            params.doChainRef = true;
            params.sigmaHSRef = 6;
//            params.BDtol = 1e-11;
        }
        // Import Params
        final ChemForm[] chemForm = params.chemForm;
        final int nPoints = params.nPoints;
        int[] t = params.types;
        if (t == null) {
            if (chemForm.length > 1) {
                throw new RuntimeException("If you have a mixture, you must specify types");
            }
            t = new int[nPoints];
        }
        else if (t.length != nPoints) {
            throw new RuntimeException("types must be of length nPoints");
        }

        final int[] nTypes = new int[chemForm.length];
        final int[] types = t;
        for (int i=0; i<t.length; i++) {
            nTypes[t[i]]++;
        }
        final int nDer = params.nDer;
        final double temperatureK = params.temperature;
        final long steps = params.numSteps;

        double refFrac = params.refFrac;
        double sigmaHSRef = params.sigmaHSRef;
        int[] seed = params.seed;

        boolean dorefpref = params.dorefpref;
        boolean doChainRef = params.doChainRef;

        double BDtol = params.BDtol;


        // Set Number of Blocks
        final long numBlocks = 1000;

        // Set Big Decimal Acceptance Fraction
        final double BDAccFrac = 0.1;

        // Check Params
        if (chemForm.length != nTypes.length) throw new RuntimeException("chemFrom and nTypes lengths are unequal!");

        if (chemForm.length > 1 && Arrays.stream(nTypes).sum() != nPoints)
            throw new RuntimeException("nPoints and nTypes do not match!");

        boolean isMixture = (nTypes.length > 1);

        // Check if Pure or Mixture
        if (isMixture) {
            for (int i = 0; i < nTypes.length; i++) {
                if (nTypes[i] == nPoints) isMixture = false;
            }
        }

        // Set Simulation Temperature
        double temperature = Kelvin.UNIT.toSim(temperatureK);

        // Evaluate Hard Sphere Coefficient
        double vhs = (4.0 / 3.0) * Math.PI * sigmaHSRef * sigmaHSRef * sigmaHSRef;
//        vhs = 2 * sigmaHSRef; //1-D vhs
        final double HSBn = doChainRef ? SpecialFunctions.factorial(nPoints) / 2 * Math.pow(vhs, nPoints - 1) : Standard.BHS(nPoints, sigmaHSRef);
        System.out.println("Chemform Length:" + chemForm.length);
        System.out.println(Kelvin.UNIT.toSim(98));
        // Print Pretext
        if (!isMixture) {
            ChemForm chemFormPure = chemForm[0];
            if (nTypes.length > 1) {
                for (int i = 0; i < nTypes.length; i++) {
                    if (nTypes[i] == nPoints) chemFormPure = chemForm[i];
                }
            }

            System.out.println("Overlap sampling for TraPPE " + chemFormPure + " at " + temperatureK + " K " + "for B" + nPoints + " and " + nDer + " derivatives");
        } else {
            String nTstr = "{";
            for (int i = 0; i < nTypes.length; i++) {
                if (nTypes[i] != 0) nTstr += ((nTstr == "{") ? "" : ",") + nTypes[i];
            }
            nTstr += "}";

            String CFstr = "";
            for (int i = 0; i < chemForm.length; i++) {
                if (nTypes[i] != 0) CFstr += chemForm[i] + " ";
            }

            System.out.println("Overlap sampling for TraPPE " + CFstr + " " + nTstr + " Mixture at " + temperatureK + " K " + "for B" + nPoints + " and " + nDer + " derivatives");
        }

        System.out.println("Reference diagram: B" + nPoints + " for hard spheres with diameter " + sigmaHSRef + " Angstroms");

        System.out.println("  B" + nPoints + "HS: " + HSBn);

        // Set up Space
        Space space = Space3D.getInstance();

        // Setting up Reference Cluster Mayer Function
        MayerFunction fRefPos = new MayerFunction() {
            public void setBox(Box box) {
            }

            public double f(IMoleculeList pair, double r2, double beta) {
                return r2 < sigmaHSRef * sigmaHSRef ? 1 : 0;
            }
        };

        // Setting up Reference Cluster
        MayerHardSphere fRef = new MayerHardSphere(sigmaHSRef);
        ClusterAbstract refCluster = doChainRef ? new ClusterChainHS(nPoints, fRefPos) : new ClusterWheatleyHS(nPoints, fRef);
        refCluster.setTemperature(temperature);

        // Setting up Target Cluster Mayer Function
        ClusterAbstractMultivalue targetClusterRigid = null;

        ClusterWheatleySoftDerivativesBD targetClusterBDRigid = null;
        SpeciesManager.Builder sb = SpeciesManager.builder();

        boolean anyPolar = false;
        TraPPEParams[] TPList = new TraPPEParams[chemForm.length];
        for (int i = 0; i < TPList.length; i++) {
            TPList[i] = new TraPPEParams(space, chemForm[i]);
        }

        for (int i = 0; i < chemForm.length; i++) {
            sb.addSpecies(TPList[i].species);
        }
        SpeciesManager sm = sb.build();

        IPotentialMolecular[][] potentials = new IPotentialMolecular[nTypes.length][nTypes.length];
        for (int i = 0; i < chemForm.length; i++) {
            TraPPEParams TPi = TPList[i];
            TPi.buildPotentials(sm);
            PotentialMoleculePair PGii = TPi.potentialGroup;

            P2PotentialGroupBuilder.ModelParams MPi = new P2PotentialGroupBuilder.ModelParams(TPi.atomTypes, TPi.sigma, TPi.epsilon, TPi.charge);
            potentials[i][i] = PGii;

            anyPolar = (anyPolar || TPi.polar);

            for (int j = i + 1; j < chemForm.length; j++) {

                TraPPEParams TPj = TPList[j];

                P2PotentialGroupBuilder.ModelParams MPj = new P2PotentialGroupBuilder.ModelParams(TPj.atomTypes, TPj.sigma, TPj.epsilon, TPj.charge);

                PotentialMoleculePair PGij = P2PotentialGroupBuilder.P2PotentialGroupBuilder(space, sm, MPi, MPj);

                potentials[i][j] = potentials[j][i] = PGij;

            }
        }

        MayerMix fTarget = new MayerMix(potentials);

        //flex moves
        PotentialMasterBonding.FullBondingInfo bondingInfo = new PotentialMasterBonding.FullBondingInfo(sm);
        int nSpheres = TPList[0].species.getAtomTypes().size();
        boolean isFlex = TPList[0].isFlex && (params.diagram == null|| !params.diagram.equals("BC"));

        System.out.println("isFlex = " + isFlex);
        System.out.println("Diagram " + params.diagram);
        VirialDiagrams Diagrams = new VirialDiagrams(nPoints, false, isFlex);

        Diagrams.setDoReeHoover(false);
        ClusterAbstract targetCluster = Diagrams.makeVirialCluster(fTarget);
        ClusterSumShell[] targetDiagrams = new ClusterSumShell[0];
        int[] targetDiagramNumbers = new int[0];
        boolean[] diagramFlexCorrection = null;

        if (TPList[0].isFlex) {
            Set<Graph> singleGraphs = Diagrams.getMSMCGraphs(true, false);
            Map<Graph,Graph> cancelMap = Diagrams.getCancelMap();
            if(params.diagram != null && !params.diagram.equals("BC")) {
                int iGraph = 0;
                for (Graph g : singleGraphs) {
                    if(params.diagram.equals(g.getStore().toNumberString() + "c")) {
                        targetCluster = Diagrams.makeVirialCluster(g, fTarget);
                        System.out.print(iGraph+" ("+g.coefficient()+") "+g.getStore().toNumberString()); // toNumberString: its corresponding number
                        Graph cancelGraph = cancelMap.get(g);
                        Set<Graph> gSplit = Diagrams.getSplitDisconnectedVirialGraphs(cancelGraph);

                        System.out.print(" - "+VirialAlkane.getSplitGraphString(gSplit, Diagrams, false));

                    }
                    iGraph++;
                }
            }
            else {
                targetDiagramNumbers = new int[targetDiagrams.length];
                System.out.println("individual clusters:");
                int iGraph = 0;
                diagramFlexCorrection = new boolean[targetDiagrams.length];
                for (Graph g : singleGraphs) {
                    System.out.print(iGraph + " (" + g.coefficient() + ") " + g.getStore().toNumberString()); // toNumberString: its corresponding number
//                    targetDiagramNumbers[iGraph] = Integer.parseInt(g.getStore().toNumberString());

                    Graph cancelGraph = cancelMap.get(g);
                    if (cancelGraph != null) {
                        diagramFlexCorrection[iGraph] = true;
                        Set<Graph> gSplit = Diagrams.getSplitDisconnectedVirialGraphs(cancelGraph);

                        System.out.print(" - " + VirialAlkane.getSplitGraphString(gSplit, Diagrams, false));

                    }
                    System.out.println();
                    iGraph++;
                }
            }
            System.out.println();
            Set<Graph> disconnectedGraphs = Diagrams.getExtraDisconnectedVirialGraphs();
            if (disconnectedGraphs.size() > 0) {
                System.out.println("extra clusters:");

                for (Graph g : disconnectedGraphs) {
                    Set<Graph> gSplit = Diagrams.getSplitDisconnectedVirialGraphs(g);
                    System.out.println(g.coefficient()+" "+VirialAlkane.getSplitGraphString(gSplit, Diagrams, true));
                }
                System.out.println();
            }
        }

        targetCluster.setTemperature(temperature);
        for (int i=0; i<targetDiagrams.length; i++) {
            targetDiagrams[i].setTemperature(temperature);
        }


        //P3 bond angle
        for (TraPPEParams TP : TPList) {
            if (TP.theta_eq != null) {
                P3BondAngle[] p3 = new P3BondAngle[TP.theta_eq.length]; //declaration, instatation
                for (int j = 0; j < TP.theta_eq.length; j++) {
                    p3[j] = new P3BondAngle(TP.theta_eq[j], TP.k_theta[j]);

                    bondingInfo.setBondingPotentialTriplet(TP.species, p3[j], Collections.singletonList(TP.triplets[j]));

                }

            }
            //dihedral stuff
            if (TP.a != null) {
                if (TP.a.length == 4){
                    P4BondTorsion[] p4 = new P4BondTorsion[TP.a.length];

                    for (int i=0; i < TP.a.length; i++) {
                        p4[i] = new P4BondTorsion(space, TP.a[i][0], TP.a[i][1], TP.a[i][2], TP.a[i][3]);
                        bondingInfo.setBondingPotentialQuad(TP.species, p4[i], Collections.singletonList(TP.quads[i]));

                    }
                }
                else{
                    P4BondTorsionAlkylBenzene[] p4 = new P4BondTorsionAlkylBenzene[TP.a.length];

                    for (int i=0; i < TP.a.length; i++) {
                        p4[i] = new P4BondTorsionAlkylBenzene(space, TP.a[i][0], TP.a[i][1]);
                        bondingInfo.setBondingPotentialQuad(TP.species, p4[i], Collections.singletonList(TP.quads[i]));

                    }
                }

    //            System.exit(0);

            }

        }


        // Setting up Target Cluster for Rigid
        MayerFunction[][] fAll = new MayerFunction[nTypes.length][nTypes.length];
        for (int i=0; i<fAll.length; i++) {
            for (int j=0; j<fAll.length; j++) {
                fAll[i][j] = fTarget;
            }
        }
        targetClusterRigid = new ClusterWheatleySoftDerivatives(nPoints, fTarget, BDtol, nDer);
        targetClusterRigid.setTemperature(temperature);

        // Setting BlockSize
        long blockSize = steps/numBlocks;
        int EqSubSteps = 1000;

        System.out.println(steps + " steps (" + numBlocks + " blocks of " + blockSize + ")");
        System.out.println("BD_Tol: " + BDtol + " BDAccFrac: " + BDAccFrac);

        // Setting up Flipping rigid, polar
        if(anyPolar && !isFlex && nPoints==2) {
            System.out.println("Performing Flipping");
            ((ClusterWheatleySoftDerivatives) targetClusterRigid).setTolerance(0);
            final int precision = -3*(int)Math.log10(BDtol);
            targetClusterBDRigid = new ClusterWheatleySoftDerivativesBD(nPoints,fTarget,precision,nDer);
            targetClusterBDRigid.setTemperature(temperature);
            ((ClusterWheatleySoftDerivatives) targetClusterRigid).setDoCaching(false);
            targetClusterBDRigid.setDoCaching(false);
            targetClusterRigid = new ClusterCoupledFlippedMultivalue(targetClusterRigid, targetClusterBDRigid, space, 3, nDer, BDtol);
        }
        //flipping for flexible polar B2
        else if(isFlex && nPoints==2 && anyPolar){
            ((ClusterSum)targetCluster).setCaching(false);

            targetCluster = new ClusterCoupledFlipped(targetCluster, space, 3);

        }
        else if(anyPolar && isFlex && nPoints > 2 && params.diagram != null && !params.diagram.equals("BC") ){
            int[][] flipPoints = Diagrams.getFlipPointsforDiagram(params.diagram);
            ((ClusterSum)targetCluster).setCaching(false);

            targetCluster = new ClusterCoupledFlippedPoints(targetCluster, space, flipPoints, 3);

        }


        // Setting up Simulation
        IPotential2[][] potential = new IPotential2[sm.getAtomTypeCount()][sm.getAtomTypeCount()];
        boolean needIntraPotentials = false;

        int[] myNTypes = nTypes.clone();
        if (isFlex) {
            myNTypes[types[0]]++;
        }
        SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space, sm, myNTypes, temperature, refCluster, isFlex ? targetCluster : targetClusterRigid);
        if(seed!=null)sim.setRandom(new RandomMersenneTwister(seed));
        System.out.println("random seeds: "+ Arrays.toString(seed==null?sim.getRandomSeeds():seed));

        if(!isFlex) {
            for (TraPPEParams TP : TPList){
                if(TP.isFlex){
                    sim.setBondingInfo(bondingInfo);
                    sim.setIntraPairPotentials(potential);
                    needIntraPotentials = true;
                    break;
                }

            }
            if (targetClusterRigid instanceof ClusterCoupledFlippedMultivalue) {
                ((ClusterCoupledFlippedMultivalue) targetClusterRigid).setBDAccFrac(BDAccFrac, sim.getRandom());
            } else {
                ((ClusterWheatleySoftDerivatives) targetClusterRigid).setBDAccFrac(BDAccFrac, sim.getRandom());
                ((ClusterWheatleySoftDerivatives) targetClusterRigid).setNumBDCheckBins(8);
            }
            // Adding derivative clusters to simulation
            ClusterMultiToSingle[] primes = new ClusterMultiToSingle[nDer];
            for(int m=0;m<primes.length;m++){
                primes[m] = new ClusterMultiToSingle(targetClusterRigid, m + 1);
            }
            sim.setExtraTargetClusters(primes);

        }
        else{
            sim.setExtraTargetClusters(targetDiagrams);
//            sim.setDoWiggle(nSpheres > 2);

            sim.setBondingInfo(bondingInfo);
            needIntraPotentials = true;
            sim.setIntraPairPotentials(potential);
        }
        if (needIntraPotentials) {
            for (TraPPEParams TP : TPList){
                if (!TP.isFlex){
                    continue;
                }
                for (AtomType T1 : TP.species.getAtomTypes()){
                    for (AtomType T2 : TP.species.getAtomTypes()){
                        potential[T1.getIndex()][T2.getIndex()] = TP.potentialGroup.getAtomPotentials()[T1.getIndex()][T2.getIndex()];
                    }
                }

            }

        }
        // Initialize Simulation
        sim.init();

        if (isFlex && isMixture) {
            // we always need to do this because (if nothing else) our
            // alternate root will be out of order
            sim.box[0].setTypes(types);
            sim.box[1].setTypes(types);
        }

        // Set Position Definitions
        sim.box[0].setPositionDefinition(new MoleculePositionCOM(space));
        sim.box[1].setPositionDefinition(new MoleculePositionCOM(space));
        int[] constraintMap = new int[nPoints+1];
        MCMoveClusterAngleGeneral mcMoveAngle = null;
        MCMoveClusterAngleGeneral mcMoveAngle1 = null;
        MCMoveClusterAngleGeneral mcMoveAngle_oneSide = null;
        MCMoveClusterAngleGeneral mcMoveAngle1_oneSide = null;

        if (isFlex) {
            for (int i=0; i<nPoints; i++) {
                constraintMap[i] = i;
            }
            constraintMap[nPoints] = 0;
            ((MCMoveClusterMoleculeMulti)sim.mcMoveTranslate[1]).setConstraintMap(constraintMap);
            ((MCMoveClusterRotateMoleculeMulti)sim.mcMoveRotate[0]).setConstraintMap(constraintMap);
            ((MCMoveClusterRotateMoleculeMulti)sim.mcMoveRotate[1]).setConstraintMap(constraintMap);


        }
        for (TraPPEParams TP : TPList) {

            if(TP.theta_eq != null) {
                mcMoveAngle = new MCMoveClusterAngleGeneral(sim.integrators[0].getPotentialCompute(), space,TP.species, TP.bonding, false, TP.triplets, sim.getRandom(), 0.1);
                mcMoveAngle.setStepSizeMax(0.6);
                mcMoveAngle.setConstraintMap(constraintMap);
                sim.integrators[0].getMoveManager().addMCMove(mcMoveAngle);
                mcMoveAngle1 = new MCMoveClusterAngleGeneral(sim.integrators[1].getPotentialCompute(), space,TP.species, TP.bonding, false, TP.triplets, sim.getRandom(), 0.1);
                mcMoveAngle1.setConstraintMap(constraintMap);
                sim.integrators[1].getMoveManager().addMCMove(mcMoveAngle1);
                mcMoveAngle1.setStepSizeMax(0.6);
                ((MCMoveStepTracker)mcMoveAngle1.getTracker()).setNoisyAdjustment(true);

                if (TP.bonding.length > 3) {
                    mcMoveAngle_oneSide = new MCMoveClusterAngleGeneral(sim.integrators[0].getPotentialCompute(), space, TP.species, TP.bonding, true, TP.triplets, sim.getRandom(), 0.1);
                    mcMoveAngle_oneSide.setStepSizeMax(0.6);
                    mcMoveAngle_oneSide.setConstraintMap(constraintMap);
                    sim.integrators[0].getMoveManager().addMCMove(mcMoveAngle_oneSide);
                    mcMoveAngle1_oneSide = new MCMoveClusterAngleGeneral(sim.integrators[1].getPotentialCompute(), space, TP.species, TP.bonding, true, TP.triplets, sim.getRandom(), 0.1);
                    sim.integrators[1].getMoveManager().addMCMove(mcMoveAngle1_oneSide);
                    mcMoveAngle1_oneSide.setStepSizeMax(0.6);
                    mcMoveAngle1_oneSide.setConstraintMap(constraintMap);
                    ((MCMoveStepTracker)mcMoveAngle1_oneSide.getTracker()).setNoisyAdjustment(true);

                }


        }}
        if (doChainRef) {
            sim.integrators[0].getMoveManager().removeMCMove(sim.mcMoveTranslate[0]);

            MCMoveClusterMoleculeHSChain mcMoveHSC = new MCMoveClusterMoleculeHSChain(sim.getRandom(), sim.box[0], sigmaHSRef);
            if(isFlex) {
                mcMoveHSC.setConstraintMap(constraintMap);
            }
            sim.integrators[0].getMoveManager().addMCMove(mcMoveHSC);
            sim.accumulators[0].setBlockSize(1);
        }
        ((MCMoveStepTracker)sim.mcMoveTranslate[1].getTracker()).setNoisyAdjustment(true);
        ((MCMoveStepTracker)sim.mcMoveRotate[1].getTracker()).setNoisyAdjustment(true);

        // create the intramolecular potential here, add to it and add it to
        // the potential master if needed
        MCMoveClusterTorsion mcMoveTorsion = null;
        MCMoveClusterTorsion mcMoveTorsion1 = null;
        for (TraPPEParams TP : TPList) {
            if (TP.a != null) {
                mcMoveTorsion = new MCMoveClusterTorsion(sim.integrators[0].getPotentialCompute(), space,TP.species, TP.bonding, TP.quads, sim.getRandom(), 1);
                mcMoveTorsion.setStepSizeMax(2);
                mcMoveTorsion.setConstraintMap(constraintMap);
                sim.integrators[0].getMoveManager().addMCMove(mcMoveTorsion);
                mcMoveTorsion1 = new MCMoveClusterTorsion(sim.integrators[1].getPotentialCompute(), space,TP.species, TP.bonding, TP.quads, sim.getRandom(), 1);
                mcMoveTorsion1.setConstraintMap(constraintMap);
                ((MCMoveStepTracker)mcMoveTorsion1.getTracker()).setNoisyAdjustment(true);

                sim.integrators[1].getMoveManager().addMCMove(mcMoveTorsion1);
                mcMoveTorsion1.setStepSizeMax(2);
            }
        }
        if (isFlex) {
            double minValue = 1e-10;
            ((ClusterWeightAbs)sim.getSampleClusters()[1]).setMinValue(minValue);

            for (int i = 0; i < 100000 && sim.box[1].getSampleCluster().value(sim.box[1]) <= minValue; i++) {
                sim.integrators[1].doStep();
            }
            if (sim.box[1].getSampleCluster().value(sim.box[1]) <= minValue) {
                throw new RuntimeException("Could not find appropriate starting config");
            }
            System.out.println("It took " + sim.integrators[1].getStepCount() + " steps to find a starting config.");
            ((ClusterWeightAbs) sim.getSampleClusters()[1]).setMinValue(0);
        }
//        System.out.println("energy: "+sim.integrators[0].getPotentialCompute().computeAll(false));
//        for (int i = 0; i < sim.box[0].getLeafList().size(); i++){
//            System.out.println("position of " + i + ": " + sim.box[0].getLeafList().get(i).getPosition());
//        }

        // Run with Graphics
        if (false) {
            sim.box[0].getBoundary().setBoxSize(space.makeVector(new double[]{10,10,10}));
            sim.box[1].getBoundary().setBoxSize(space.makeVector(new double[]{10,10,10}));
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
            DisplayBox displayBox0 = simGraphic.getDisplayBox(sim.box[0]);
            DisplayBox displayBox1 = simGraphic.getDisplayBox(sim.box[1]);
//            displayBox0.setPixelUnit(new Pixel(300.0/size));
//            displayBox1.setPixelUnit(new Pixel(300.0/size));
            displayBox0.setShowBoundary(false);
            displayBox1.setShowBoundary(false);
            ((DisplayBoxCanvasG3DSys)displayBox0.canvas).setBackgroundColor(Color.WHITE);
            ((DisplayBoxCanvasG3DSys)displayBox1.canvas).setBackgroundColor(Color.WHITE);


            ColorSchemeRandomByMolecule colorScheme = new ColorSchemeRandomByMolecule(sim.getSpeciesManager(), sim.box[0], sim.getRandom());
            displayBox0.setColorScheme(colorScheme);
            colorScheme = new ColorSchemeRandomByMolecule(sim.getSpeciesManager(), sim.box[1], sim.getRandom());
            displayBox1.setColorScheme(colorScheme);
            simGraphic.makeAndDisplayFrame();

            sim.integratorOS.setNumSubSteps(1000);
            sim.setAccumulatorBlockSize(1000);

            // if running interactively, set filename to null so that it doens't read
            // (or write) to a refpref file
            sim.initRefPref(null, 100, false);
            sim.equilibrate(null, 200, false);
            sim.getController().addActivity(new ActivityIntegrate(sim.integratorOS));
            return;
        }

        // Setting up Equilibration
        sim.integratorOS.setNumSubSteps(EqSubSteps);
        sim.integratorOS.setAggressiveAdjustStepFraction(true);

        // Start timing
        long t1 = System.currentTimeMillis();

        if (refFrac >= 0) {
            sim.integratorOS.setRefStepFraction(refFrac);
            sim.integratorOS.setAdjustStepFraction(false);
        }

        System.out.println();
        String refFileName = null;
        System.out.println("iscommandline:"+isCommandline);
        if (isCommandline) {
            // if running interactively, don't use the file
            String tempString = ""+temperatureK;
            if (temperatureK == (int)temperatureK) {
                // temperature is an integer, use "200" instead of "200.0"
                tempString = ""+(int)temperatureK;
            }
            refFileName = "refpref_"+"_"+nPoints+"_"+tempString+"K"+ Arrays.toString(params.chemForm);
            if(params.diagram != null){
                refFileName += "_"+params.diagram;
            }
        }
        sim.integrators[1].getEventManager().addListener(new IntegratorListener() {
            @Override
            public void integratorStepFinished(IntegratorEvent e) {
                if (sim.integrators[1].getStepCount() % 10000 != 0) return;
                System.out.println(sim.integrators[1].getStepCount());
                PDBWriter writer = new PDBWriter();
                writer.setBox(sim.box[1]);
                writer.setFileName(sim.integrators[1].getStepCount() + ".pdb");
                writer.actionPerformed();
            }
        });
        // Equilibrate
        sim.initRefPref(refFileName, (steps / EqSubSteps) / 20);
        sim.equilibrate(refFileName, (steps / EqSubSteps) / 10);

        System.out.println("equilibration finished");

        if(dorefpref){
            long t2 = System.currentTimeMillis();
            System.out.println("time: "+(t2-t1)/1000.0);
            return;
        }


        // Setting up Production Run
        sim.integratorOS.setNumSubSteps((int) blockSize);
        sim.setAccumulatorBlockSize(blockSize);

        if (doChainRef) sim.accumulators[0].setBlockSize(1);
        for (int i = 0; i < 2; i++) {
            if (i > 0 || !doChainRef) System.out.println("MC Move step sizes " + sim.mcMoveTranslate[i].getStepSize());
        }

        // Production Run
        ActivityIntegrate ai = new ActivityIntegrate(sim.integratorOS, numBlocks);
        sim.getController().runActivityBlocking(ai);

        System.out.println();


        // Print BD and Flip Stats
        if(!isFlex) {
            // Print Simulation Output
            System.out.println("final reference step fraction " + sim.integratorOS.getIdealRefStepFraction());
            System.out.println("actual reference step fraction " + sim.integratorOS.getRefStepFraction());

            String[] extraNames = new String[nDer];
            for (int i = 1; i <= nDer; i++) {
                extraNames[i - 1] = "derivative " + i;
            }
            sim.printResults(HSBn, extraNames);

            if (targetClusterRigid instanceof ClusterWheatleySoftDerivatives) {
                System.out.println("SoftBDcount: " + ((ClusterWheatleySoftDerivatives) targetClusterRigid).getSoftBDcount() + " SoftBDfrac: " + ((ClusterWheatleySoftDerivatives) targetClusterRigid).getSoftBDfrac() + " Softcount: " + ((ClusterWheatleySoftDerivatives) targetClusterRigid).getSoftcount());
            } else if (targetClusterRigid instanceof ClusterCoupledFlippedMultivalue) {
                ClusterCoupledFlippedMultivalue foo = (ClusterCoupledFlippedMultivalue) targetClusterRigid;
                System.out.println("BDcount: " + foo.getBDcount() + " BDfrac: " + foo.getBDfrac() + " totBDcount: " + foo.getBDtotcount());
                System.out.println("FlipCount: " + foo.getflipcount() + " Flipfrac: " + foo.getflipfrac() + " FlipTotcount: " + foo.gettotcount());
            }
        }
        long t2 = System.currentTimeMillis();
        System.out.println("time: " + (t2 - t1) / 1000.0);
        if(isFlex){
            System.out.println("Angle move acceptance "+ mcMoveAngle1.getTracker().acceptanceRatio() + " " + mcMoveAngle.getTracker().acceptanceRatio());
            if (nSpheres > 3) {
                for (int i=0; i<TPList[0].a.length; i++) {
                System.out.println("Torsion move acceptance "+mcMoveTorsion.getTracker().acceptanceRatio());}
            }

            System.out.println("final reference step frequency "+sim.integratorOS.getIdealRefStepFraction());
            System.out.println("actual reference step frequency "+sim.integratorOS.getRefStepFraction());
            String[] extraNames = new String[targetDiagrams.length];
            for (int i=0; i<targetDiagrams.length; i++) {
                String n = "";
                if (targetDiagramNumbers[i] < 0) {
                    n = "diagram " + (-targetDiagramNumbers[i]) + "bc";
                } else {
                    n = "diagram " + targetDiagramNumbers[i];

                    if (diagramFlexCorrection[i]) {
                        n += "c";
                    }
                }
                extraNames[i] = n;
            }
            sim.printResults(HSBn, extraNames);
        }

    }

    enum ChemForm {
        N2, O2, CO2, NH3, CH4, CH3OH, ethanol, propan1ol, propan2ol, isobutanol, benzene, ethane, propane, butane, methane, ethene, propene, butadiene, toluene,
        ethylbenzene, oxylene, pxylene, mxylene, water
    }

    /**
     * Inner class for parameters
     */
    public static class VirialParam extends ParameterBase {
        // don't change these
        public ChemForm[] chemForm = {ChemForm.N2};
        public int nPoints = 2;
        public int[] types = null;  // null works for pure component
        public int nDer = 3;
        public double temperature = 400;
        public long numSteps = 1000000;

        public double refFrac = -1;
        public double sigmaHSRef = 5;
        public int[] seed = null;

        public boolean dorefpref = true;
        public boolean doChainRef = true;
        public String diagram = null;
        public double BDtol = 1e-11;

    }

    /**
     * Inner class for TraPPE Parameters
     */
    public static class TraPPEParams{

        public final AtomType[] atomTypes;
        public final double[] sigma;
        public final double[] epsilon;
        public final double[] charge;
        public double[] k_theta;
        public double[] theta_eq;
        public double[][] a;
        public int[][] triplets;
        public int[][] quads;
        public IntArrayList[] bonding;
        public final ISpecies species;
        public PotentialMoleculePair potentialGroup;
        protected static Element elementM = new ElementSimple("M", 0.0);
        protected boolean isFlex;
        protected boolean polar;
        protected final ChemForm chemForm;
        protected final Space space;

        //Set up computing the boolean. It is hard coded for now.

        public TraPPEParams(Space space, ChemForm chemForm){
            this.chemForm = chemForm;
            this.space = space;

            if(chemForm == ChemForm.N2) {

                //Atoms in Compound
                AtomType typeM = new AtomType(elementM);
                AtomType typeN = new AtomType(Nitrogen.INSTANCE);

                atomTypes = new AtomType[]{typeM,typeN};

                //TraPPE Parameters
                double bondLength = 1.10; // Angstrom
                double sigmaN = 3.31; // Angstrom
                double epsilonN = Kelvin.UNIT.toSim(36.0);
                double qN = Electron.UNIT.toSim(-0.482);
                double sigmaM = 0.0; // Angstrom
                double epsilonM = Kelvin.UNIT.toSim(0.0);
                double qM = Electron.UNIT.toSim(0.964);
                isFlex = false; //rigid
                //Construct Arrays
                sigma = new double[]{sigmaM, sigmaN};
                epsilon = new double[]{epsilonM, epsilonN};
                charge = new double[]{qM, qN};

                //Get Coordinates
                Vector3D posM = new Vector3D(new double[]{0, 0, 0});
                Vector3D posN1 = new Vector3D(new double[]{-bondLength / 2, 0, 0});
                Vector3D posN2 = new Vector3D(new double[]{+bondLength / 2, 0, 0});

                //Set Geometry
                species = new SpeciesBuilder(space)
                        .addAtom(typeM, posM, "M")
                        .addAtom(typeN, posN1, "N1")
                        .addAtom(typeN, posN2, "N2")
                        .build();
            }

            else if (chemForm == ChemForm.O2) {

                //Atoms in Compound
                AtomType typeM = new AtomType(elementM);
                AtomType typeO = new AtomType(Oxygen.INSTANCE);

                atomTypes = new AtomType[]{typeM,typeO};
                isFlex = false;
                //TraPPE Parameters
                double bondLength = 1.210; // Angstrom
                double sigmaO = 3.020; // Angstrom
                double epsilonO = Kelvin.UNIT.toSim(49.0);
                double qO = Electron.UNIT.toSim(-0.113);
                double sigmaM = 0.0; // Angstrom
                double epsilonM = Kelvin.UNIT.toSim(0.0);
                double qM = Electron.UNIT.toSim(0.226);

                //Construct Arrays
                sigma = new double[]{sigmaM, sigmaO};
                epsilon = new double[]{epsilonM, epsilonO};
                charge = new double[]{qM, qO};

                //Get Coordinates
                Vector3D posM = new Vector3D(new double[]{0, 0, 0});
                Vector3D posO1 = new Vector3D(new double[]{-bondLength / 2, 0, 0});
                Vector3D posO2 = new Vector3D(new double[]{+bondLength / 2, 0, 0});

                //Set Geometry
                species = new SpeciesBuilder(space)
                        .addAtom(typeM, posM, "M")
                        .addAtom(typeO, posO1, "O1")
                        .addAtom(typeO, posO2, "O2")
                        .build();
            }

            else if (chemForm == ChemForm.CO2) {

                //Atoms in Compound
                AtomType typeC = new AtomType(Carbon.INSTANCE);
                AtomType typeO = new AtomType(Oxygen.INSTANCE);

                atomTypes = new AtomType[]{typeC,typeO};
                isFlex = false;
                //TraPPE Parameters
                double bondLengthCO = 1.160; // Angstrom
                double sigmaC = 2.800; // Angstrom
                double epsilonC = Kelvin.UNIT.toSim(27.0);
                double qC = Electron.UNIT.toSim(0.700);
                double sigmaO = 3.050; // Angstrom
                double epsilonO = Kelvin.UNIT.toSim(79.0);
                double qO = Electron.UNIT.toSim(-0.350);

                //Construct Arrays
                sigma = new double[]{sigmaC, sigmaO};
                epsilon = new double[]{epsilonC, epsilonO};
                charge = new double[]{qC, qO};

                //Get Coordinates
                Vector3D posC = new Vector3D(new double[]{0, 0, 0});
                Vector3D posO1 = new Vector3D(new double[]{-bondLengthCO, 0, 0});
                Vector3D posO2 = new Vector3D(new double[]{+bondLengthCO, 0, 0});

                //Set Geometry
                species = new SpeciesBuilder(space)
                        .addAtom(typeC, posC, "C")
                        .addAtom(typeO, posO1, "O1")
                        .addAtom(typeO, posO2, "O2")
                        .build();
            }
            else if (chemForm == ChemForm.NH3) {

                //Atom in Compound
                AtomType typeN = new AtomType(Nitrogen.INSTANCE);
                AtomType typeH = new AtomType(Hydrogen.INSTANCE);
                AtomType typeM = new AtomType(elementM);

                atomTypes = new AtomType[]{typeN,typeH,typeM};

                polar = true;
                isFlex = false;
                //TraPPE Parameters
                double bondLengthNH = 1.012; // Angstrom
                double bondLengthNM = 0.080; // Angstrom
                double thetaHNM = Degree.UNIT.toSim(67.9) ;
                double thetaHNH = Degree.UNIT.toSim(106.7);
                double thetaHNHxy = Degree.UNIT.toSim(60);
                double sigmaN = 3.420; // Angstrom
                double epsilonN = Kelvin.UNIT.toSim(185.0);
                double qN = Electron.UNIT.toSim(0.0);
                double sigmaH = 0.0; // Angstrom
                double epsilonH = Kelvin.UNIT.toSim(0.0);
                double qH = Electron.UNIT.toSim(0.410);
                double sigmaM = 0.0; // Angstrom
                double epsilonM = Kelvin.UNIT.toSim(0.0);
                double qM = Electron.UNIT.toSim(-1.230);

                //Construct Arrays
                sigma = new double[] {sigmaN,sigmaH,sigmaM};
                epsilon = new double[] {epsilonN,epsilonH,epsilonM};
                charge = new double[]{qN, qH, qM};

                //Get Coordinates
                Vector3D posN = new Vector3D(new double[]{0, 0, 0});
                Vector3D posH1 = new Vector3D(new double[]{bondLengthNH * Math.sin(thetaHNM), 0, -bondLengthNH * Math.cos(thetaHNM)});
                Vector3D posH2 = new Vector3D(new double[]{-bondLengthNH * Math.sin(thetaHNM) * Math.cos(thetaHNHxy), bondLengthNH * Math.sin(thetaHNM) * Math.sin(thetaHNHxy), -bondLengthNH * Math.cos(thetaHNM)});
                Vector3D posH3 = new Vector3D(new double[]{-bondLengthNH * Math.sin(thetaHNM) * Math.cos(thetaHNHxy), -bondLengthNH * Math.sin(thetaHNM) * Math.sin(thetaHNHxy), -bondLengthNH * Math.cos(thetaHNM)});
                Vector3D posM = new Vector3D(new double[]{0, 0, -bondLengthNM});

                //Set Geometry
                species = new SpeciesBuilder(space)
                        .addAtom(typeN, posN, "N")
                        .addAtom(typeH, posH1, "H1")
                        .addAtom(typeH, posH2, "H2")
                        .addAtom(typeH, posH3, "H3")
                        .addAtom(typeM, posM, "M")
                        .build();
            }
//            else if (chemForm == ChemForm.CH4) {
//                //TraPPE-EH
//                //Atom in Compound
//                //Avogadro 3D checked
//                AtomType typeC = new AtomType(Carbon.INSTANCE);
//                AtomType typeM = new AtomType(elementM);
//                isFlex = false;
//                atomTypes = new AtomType[]{typeC,typeM};
//
//                //TraPPE Parameters
////                double bondLengthCM = 0.55; // Angstrom
//                double sigmaC = 3.31; // Angstrom
//                double epsilonC = Kelvin.UNIT.toSim(0.01);
//                double qC = Electron.UNIT.toSim(0.0);
//                double sigmaM = 3.31; // Angstrom
//                double epsilonM = Kelvin.UNIT.toSim(15.30);
//                double qM = Electron.UNIT.toSim(0.000);
//
//                //Construct Arrays
//                sigma = new double[] {sigmaC,sigmaM};
//                epsilon = new double[] {epsilonC,epsilonM};
//                charge = new double[]{qC, qM};
//
//                //Get Coordinates
//                Vector3D posC = new Vector3D(new double[]{0, 0, 0});
//                Vector3D posM1 = new Vector3D(new double[]{0.0000000000, -0.4490790744,     -0.3175342264});
//                Vector3D posM2 = new Vector3D(new double[]{0.0000000000,      0.4490790744,     -0.3175342264});
//                Vector3D posM3 = new Vector3D(new double[]{-0.4490790744,      0.0000000000,      0.3175342264});
//                Vector3D posM4 = new Vector3D(new double[]{0.4490790744,      0.0000000000,      0.3175342264  });
//
//                //Set Geometry
//                species = new SpeciesBuilder(space)
//                        .addAtom(typeC, posC, "C")
//                        .addAtom(typeM, posM1, "M1")
//                        .addAtom(typeM, posM2, "M2")
//                        .addAtom(typeM, posM3, "M3")
//                        .addAtom(typeM, posM4, "M4")
//                        .build();
//            }
            else if (chemForm == ChemForm.CH3OH) {
                //TraPPE-UA
                //ok
                AtomType typeCH3 = new AtomType(Carbon.INSTANCE);
                AtomType typeO = new AtomType(Oxygen.INSTANCE);
                AtomType typeH = new AtomType(Hydrogen.INSTANCE);

                atomTypes = new AtomType[]{typeCH3,typeO, typeH};
                polar = true;
                isFlex = true;
                //TraPPE Parameters
                double bondLengthCH3OH = 1.43; // Angstrom
                double bondLengthOH = 0.945; //Angstrom
                double theta_CH3OH = Degree.UNIT.toSim(108.5) ;
                double sigmaCH3 = 3.75; // Angstrom
                double epsilonCH3 = Kelvin.UNIT.toSim(98);
                double qCH3 = Electron.UNIT.toSim(0.265);
                double sigmaO = 3.02; // Angstrom
                double epsilonO = Kelvin.UNIT.toSim(93);
                double qO = Electron.UNIT.toSim(-0.700);
                double sigmaH = 0.0; // Angstrom
                double epsilonH = Kelvin.UNIT.toSim(0.0);
                double qH = Electron.UNIT.toSim(0.435);
                double k_thetaCH3OH = Kelvin.UNIT.toSim(55400.0);

                double xCH3 = -bondLengthCH3OH, xH = -bondLengthOH*Math.cos(theta_CH3OH);
                double yH = bondLengthOH*Math.sin(theta_CH3OH);
                double mm = typeCH3.getMass() + typeO.getMass() + typeH.getMass();
                double xCOM = (xCH3*typeCH3.getMass() + xH * typeH.getMass()) / mm;
                double yCOM = yH*typeH.getMass() / mm;

                //Construct Arrays
                sigma = new double[] {sigmaCH3,sigmaO, sigmaH};
                epsilon = new double[] {epsilonCH3,epsilonO, epsilonH};
                charge = new double[]{qCH3, qO, qH};
                theta_eq = new double[]{theta_CH3OH};
                k_theta = new double[]{k_thetaCH3OH};
                triplets = new int[][]{{0, 1, 2}};
                bonding = new IntArrayList[3];
                bonding[0] = new IntArrayList(new int[]{1});
                bonding[1] = new IntArrayList(new int[]{0,2});
                bonding[2] = new IntArrayList(new int[]{1});

                //Get Coordinates
                Vector posCH3 = Vector.of(xCH3-xCOM, -yCOM, 0);
                Vector posO = Vector.of(-xCOM, -yCOM, 0);
                Vector posH = Vector.of(xH-xCOM, yH-yCOM, 0);

                //Set Geometry
                species = new SpeciesBuilder(space)
                        .addAtom(typeCH3, posCH3, "CH3")
                        .addAtom(typeO, posO, "O")
                        .addAtom(typeH, posH, "H")
                        .build();
            }
            else if (chemForm == ChemForm.ethanol) {
                //TraPPE-UA
                //checked
                //Atom in Compound
                AtomType typeCH3 = new AtomType(Carbon.INSTANCE);
                AtomType typeCH2 = new AtomType(Carbon.INSTANCE);
                AtomType typeO = new AtomType(Oxygen.INSTANCE);
                AtomType typeH = new AtomType(Hydrogen.INSTANCE);

                atomTypes = new AtomType[]{typeCH3, typeCH2, typeO, typeH};
                polar = true;
                isFlex = true;
                //TraPPE Parameters
                double bondLengthCC= 1.54; //Angstrom
                double bondLengthCH3OH = 1.43; // Angstrom
                double bondLengthOH = 0.945; //Angstrom
                double theta_CH3OH = Degree.UNIT.toSim(108.5) ;
                double theta_CCOH = Degree.UNIT.toSim(109.47) ;

                double sigmaCH3 = 3.75; // Angstrom
                double epsilonCH3 = Kelvin.UNIT.toSim(98);
                double qCH3 = Electron.UNIT.toSim(0.0);
                double sigmaCH2 = 3.95; // Angstrom
                double epsilonCH2 = Kelvin.UNIT.toSim(46);
                double qCH2 = Electron.UNIT.toSim(0.265);

                double sigmaO = 3.02; // Angstrom
                double epsilonO = Kelvin.UNIT.toSim(93);
                double qO = Electron.UNIT.toSim(-0.700);
                double sigmaH = 0.0; // Angstrom
                double epsilonH = Kelvin.UNIT.toSim(0.0);
                double qH = Electron.UNIT.toSim(0.435);
                double k_thetaCH3OH = Kelvin.UNIT.toSim(55400.0);
                double k_thetaCCOH = Kelvin.UNIT.toSim(50400);
                double c0 = 0.00;
                double c1 = 209.82;
                double c2 = -29.17;
                double c3 = 187.93;

                //Construct Arrays
                sigma = new double[] {sigmaCH3,sigmaCH2,sigmaO, sigmaH};
                epsilon = new double[] {epsilonCH3,epsilonCH2, epsilonO, epsilonH};
                charge = new double[]{qCH3,qCH2, qO, qH};
                theta_eq = new double[]{theta_CCOH, theta_CH3OH};
                k_theta = new double[]{k_thetaCCOH, k_thetaCH3OH};
                a = new double[][]{{c0, c1, c2, c3}};
                triplets = new int[][]{{0, 1, 2}, {1, 2, 3}};
                quads = new int[][]{{0, 1, 2, 3}};
                bonding = new IntArrayList[4];
                bonding[0] = new IntArrayList(new int[]{1});
                bonding[1] = new IntArrayList(new int[]{0,2});
                bonding[2] = new IntArrayList(new int[]{1, 3});
                bonding[3] = new IntArrayList(new int[]{2});

                double x3 = bondLengthCC - bondLengthCH3OH * Math.cos(theta_CCOH);
                double y3 = - bondLengthCH3OH * Math.sin(theta_CCOH);
                double xH = x3 + bondLengthOH * Math.cos(theta_CH3OH);
                double yH = y3 - bondLengthOH * Math.sin(theta_CH3OH);

                //Get Coordinates
                Vector3D posCH2 = new Vector3D(new double[]{0, 0,      0.0000000000 });
                Vector3D posCH3 = new Vector3D(new double[]{bondLengthCC ,     0,      0.0000000000});
                Vector3D posO = new Vector3D(new double[]{x3, y3,      0.0000000000});

                Vector3D posH = new Vector3D(new double[]{xH, yH,      0.0000000000});

                System.out.println("CH2:" + posCH2 + "CH3:" + posCH3 + "O:" + posO + "H:" + posH);

                //Set Geometry
                species = new SpeciesBuilder(space)
                        .addAtom(typeCH3, posCH3, "CH3")
                        .addAtom(typeCH2, posCH2, "CH2")
                        .addAtom(typeO, posO, "O")
                        .addAtom(typeH, posH, "H")
                        .build();
            }
            else if (chemForm == ChemForm.propan1ol) {
                //TraPPE-UA
                //Atom in Compound
                //checked
                //https://publications.lib.chalmers.se/records/fulltext/180137/local_180137.pdf
                AtomType typeCH3 = new AtomType(Carbon.INSTANCE);
                AtomType typeCH2_2 = new AtomType(Carbon.INSTANCE);
                AtomType typeCH2_3 = new AtomType(Carbon.INSTANCE);

                AtomType typeO = new AtomType(Oxygen.INSTANCE);
                AtomType typeH = new AtomType(Hydrogen.INSTANCE);

                atomTypes = new AtomType[]{typeCH3, typeCH2_2, typeCH2_3, typeO, typeH};
                polar = true;
                isFlex = true;
                //TraPPE Parameters
                double bondLengthCC= 1.54; //Angstrom
                double bondLengthCH3OH = 1.43; // Angstrom
                double bondLengthOH = 0.945; //Angstrom
                double theta_CH3OH = Degree.UNIT.toSim(108.5) ;
                double theta_CCOH = Degree.UNIT.toSim(109.47) ;
                double theta_CCC = Degree.UNIT.toSim(114) ;

                double sigmaCH3 = 3.75; // Angstrom
                double epsilonCH3 = Kelvin.UNIT.toSim(98);
                double qCH3 = Electron.UNIT.toSim(0.0);
                double sigmaCH2_2 = 3.95; // Angstrom
                double epsilonCH2_2 = Kelvin.UNIT.toSim(46);
                double qCH2_2 = Electron.UNIT.toSim(0);
                double sigmaCH2_3 = 3.95; // Angstrom
                double epsilonCH2_3 = Kelvin.UNIT.toSim(46);
                double qCH2_3 = Electron.UNIT.toSim(0.265);

                double sigmaO = 3.02; // Angstrom
                double epsilonO = Kelvin.UNIT.toSim(93);
                double qO = Electron.UNIT.toSim(-0.700);
                double sigmaH = 0.0; // Angstrom
                double epsilonH = Kelvin.UNIT.toSim(0.0);
                double qH = Electron.UNIT.toSim(0.435);
                double k_thetaCH3OH = Kelvin.UNIT.toSim(55400.0);
                double k_thetaCCOH = Kelvin.UNIT.toSim(50400);
                double k_thetaCCC = Kelvin.UNIT.toSim(62500);
                double c00 = 0.00;
                double c01 = 176.62;
                double c02 = -53.34;
                double c03 = 769.93;

                double c10 = 0.00;
                double c11 = 209.82;
                double c12 = -29.17;
                double c13 = 187.93;

                //Construct Arrays
                sigma = new double[] {sigmaCH3, sigmaCH2_2, sigmaCH2_3,sigmaO, sigmaH};
                epsilon = new double[] {epsilonCH3,epsilonCH2_2, epsilonCH2_3, epsilonO, epsilonH};
                charge = new double[]{qCH3,qCH2_2, qCH2_3, qO, qH};
                theta_eq = new double[]{theta_CCC, theta_CCOH, theta_CH3OH};
                k_theta = new double[]{k_thetaCCC, k_thetaCCOH, k_thetaCH3OH};
                a = new double[][]{{c00, c01, c02, c03}, {c10, c11, c12, c13}};
                triplets = new int[][]{{0, 1, 2}, {1, 2, 3}, {2, 3, 4}};
                quads = new int[][]{{0, 1, 2, 3}, {1, 2, 3, 4}};
                bonding = new IntArrayList[5];
                bonding[0] = new IntArrayList(new int[]{1});
                bonding[1] = new IntArrayList(new int[]{0,2});
                bonding[2] = new IntArrayList(new int[]{1,3});
                bonding[3] = new IntArrayList(new int[]{2,4});
                bonding[4] = new IntArrayList(new int[]{3});



                double x3 = bondLengthCC - bondLengthCC * Math.cos(theta_CCC);
                double y3 = bondLengthCC * Math.sin(theta_CCC);
                double xO = x3 + bondLengthCH3OH * Math.cos(theta_CCOH);
                double yO = y3 + bondLengthCH3OH * Math.sin(theta_CCOH);
                double xH = xO + bondLengthOH * Math.cos(theta_CH3OH);
                double yH = yO + bondLengthOH * Math.sin(theta_CH3OH);
                //Get Coordinates
                Vector3D posCH3 = new Vector3D(new double[]{0, 0, 0});
                Vector3D posCH2_2 = new Vector3D(new double[]{bondLengthCC ,     0,      0.0000000000});
                Vector3D posCH2_3 = new Vector3D(new double[]{x3, y3,      0.0000000000});

                Vector3D posO = new Vector3D(new double[]{xO, yO, 0});

                Vector3D posH = new Vector3D(new double[]{xH, yH, 0});
                System.out.println("CH3:" + posCH3 + "CH2_2:" + posCH2_2 + "CH2_3:" + posCH2_3 + "O:" + posO + "H:" + posH);

                //Set Geometry
                species = new SpeciesBuilder(space)
                        .addAtom(typeCH3, posCH3, "CH3")
                        .addAtom(typeCH2_2, posCH2_2, "CH2_2")
                        .addAtom(typeCH2_3, posCH2_3, "CH2_3")
                        .addAtom(typeO, posO, "O")
                        .addAtom(typeH, posH, "H")
                        .build();
            }
            else if (chemForm == ChemForm.propan2ol) {
                //TraPPE-UA
                //Atom in Compound
                //checked
                AtomType typeCH3 = new AtomType(Carbon.INSTANCE);
                AtomType typeCH = new AtomType(Carbon.INSTANCE);

                AtomType typeO = new AtomType(Oxygen.INSTANCE);
                AtomType typeH = new AtomType(Hydrogen.INSTANCE);

                atomTypes = new AtomType[]{typeCH3, typeCH, typeO, typeH};
                polar = true;
                isFlex = true;
                //TraPPE Parameters
                double bondLengthCC= 1.54; //Angstrom
                double bondLengthCH3OH = 1.43; // Angstrom
                double bondLengthOH = 0.945; //Angstrom
                double theta_CH3OH = Degree.UNIT.toSim(108.5) ;
                double theta_CCOH = Degree.UNIT.toSim(109.47) ;
                double theta_CCC = Degree.UNIT.toSim(112) ;

                double sigmaCH3 = 3.75; // Angstrom
                double epsilonCH3 = Kelvin.UNIT.toSim(98);
                double qCH3 = Electron.UNIT.toSim(0.0);
                double sigmaCH = 4.33; // Angstrom
                double epsilonCH = Kelvin.UNIT.toSim(10);
                double qCH = Electron.UNIT.toSim(0.265);

                double sigmaO = 3.02; // Angstrom
                double epsilonO = Kelvin.UNIT.toSim(93);
                double qO = Electron.UNIT.toSim(-0.700);
                double sigmaH = 0.0; // Angstrom
                double epsilonH = Kelvin.UNIT.toSim(0.0);
                double qH = Electron.UNIT.toSim(0.435);
                double k_thetaCH3OH = Kelvin.UNIT.toSim(55400.0);
                double k_thetaCCOH = Kelvin.UNIT.toSim(50400);
                double k_thetaCCC = Kelvin.UNIT.toSim(62500);
                double c00 = 215.89;
                double c01 = 197.33;
                double c02 = 31.46;
                double c03 = -173.92;


                //Construct Arrays
                sigma = new double[] {sigmaCH3, sigmaCH, sigmaO,sigmaH};
                epsilon = new double[] {epsilonCH3,epsilonCH, epsilonO, epsilonH};
                charge = new double[]{qCH3,qCH, qO, qH};
                theta_eq = new double[]{theta_CCC, theta_CCOH,  theta_CH3OH, theta_CCOH};
                k_theta = new double[]{k_thetaCCC, k_thetaCCOH, k_thetaCH3OH, k_thetaCCOH};
                a = new double[][]{{c00, c01, c02, c03}, {c00, c01, c02, c03}};
                triplets = new int[][]{{0, 1, 2}, {0, 1, 3}, {1, 3, 4}, {2, 1, 3}};
                quads = new int[][]{{0, 1, 3, 4}, {2, 1, 3, 4}};
                bonding = new IntArrayList[5];
                bonding[0] = new IntArrayList(new int[]{1});
                bonding[1] = new IntArrayList(new int[]{0,2,3});
                bonding[2] = new IntArrayList(new int[]{1});
                bonding[3] = new IntArrayList(new int[]{1,4});
                bonding[4] = new IntArrayList(new int[]{3});



                double x3 = bondLengthCC - bondLengthCC * Math.cos(theta_CCC);
                double y3 = bondLengthCC * Math.sin(theta_CCC);

                double xO = bondLengthCC - bondLengthCH3OH * Math.cos(theta_CCOH);
                double yO = - bondLengthCH3OH * Math.sin(theta_CCOH);

                double xH = xO + bondLengthOH * Math.cos(theta_CH3OH);
                double yH = yO - bondLengthOH * Math.sin(theta_CH3OH);


                //Get Coordinates
                Vector3D posCH3 = new Vector3D(new double[]{0,0,0});
                Vector3D posCH = new Vector3D(new double[]{bondLengthCC, 0, 0});
                Vector3D posCH3_3 = new Vector3D(new double[]{x3, y3, 0});

                Vector3D posO = new Vector3D(new double[]{xO, yO, 0});

                Vector3D posH = new Vector3D(new double[]{xH, yH, 0});
                System.out.println("CH3:" + posCH3 + "CH:" + posCH + "CH3_3:" + posCH3_3 + "O:" + posO + "H:" + posH);

                //Set Geometry
                species = new SpeciesBuilder(space)
                        .addAtom(typeCH3, posCH3, "CH3")
                        .addAtom(typeCH, posCH, "CH")
                        .addAtom(typeCH3, posCH3_3, "CH3_3")
                        .addAtom(typeO, posO, "O")
                        .addAtom(typeH, posH, "H")
                        .build();
            }
            else if (chemForm == ChemForm.isobutanol) {
                //TraPPE-UA
                //checked
                //Atom in Compound
                AtomType typeCH3 = new AtomType(Carbon.INSTANCE);
                AtomType typeCH = new AtomType(Carbon.INSTANCE);
                AtomType typeCH2 = new AtomType(Carbon.INSTANCE);

                AtomType typeO = new AtomType(Oxygen.INSTANCE);
                AtomType typeH = new AtomType(Hydrogen.INSTANCE);

                atomTypes = new AtomType[]{typeCH3, typeCH, typeO, typeCH2, typeH};
                polar = true;
                isFlex = true;
                //TraPPE Parameters
                double bondLengthCC= 1.54; //Angstrom
                double bondLengthCH3OH = 1.43; // Angstrom
                double bondLengthOH = 0.945; //Angstrom
                double theta_CH3OH = Degree.UNIT.toSim(108.5) ;
                double theta_CCOH = Degree.UNIT.toSim(109.47) ;
                double theta_CCC = Degree.UNIT.toSim(112) ;

                double sigmaCH3 = 3.75; // Angstrom
                double epsilonCH3 = Kelvin.UNIT.toSim(98);
                double qCH3 = Electron.UNIT.toSim(0.0);
                double sigmaCH = 4.68; // Angstrom
                double epsilonCH = Kelvin.UNIT.toSim(10);
                double qCH = Electron.UNIT.toSim(0);

                double sigmaO = 3.02; // Angstrom
                double epsilonO = Kelvin.UNIT.toSim(93);
                double qO = Electron.UNIT.toSim(-0.700);
                double sigmaH = 0.0; // Angstrom
                double epsilonH = Kelvin.UNIT.toSim(0.0);
                double qH = Electron.UNIT.toSim(0.435);
                double sigmaCH2 = 3.95; // Angstrom
                double epsilonCH2 = Kelvin.UNIT.toSim(46);
                double qCH2 = Electron.UNIT.toSim(0.265);

                double k_thetaCH3OH = Kelvin.UNIT.toSim(55400.0);
                double k_thetaCCOH = Kelvin.UNIT.toSim(50400);
                double k_thetaCCC = Kelvin.UNIT.toSim(62500);
                double c00 = 0;
                double c01 = 209.82	;
                double c02 = -29.17;
                double c03 = 187.93;
                double c10 = 0;
                double c11 = 176.62	;
                double c12 = -53.34;
                double c13 = 769.93;


                //Construct Arrays
                sigma = new double[] {sigmaCH3,sigmaCH, sigmaO, sigmaCH2, sigmaH};
                epsilon = new double[] {epsilonCH3,epsilonCH, epsilonO, epsilonCH2, epsilonH};
                charge = new double[]{qCH3,qCH, qO, qCH2, qH};
                theta_eq = new double[]{theta_CCOH, theta_CH3OH, theta_CCC, theta_CCC};
                k_theta = new double[]{k_thetaCCOH, k_thetaCH3OH, k_thetaCCC, k_thetaCCC};
                a = new double[][]{{c00, c01, c02, c03}, {c10, c11, c12, c13}, {c10, c11, c12, c13}};
                triplets = new int[][]{{4, 2, 1}, {2, 4, 5}, {2, 1, 3}, {2, 1, 0}};
                quads = new int[][]{{1, 2, 4, 5}, {0, 1, 2, 4}, {3, 1, 2, 4}};
                bonding = new IntArrayList[6];
                bonding[0] = new IntArrayList(new int[]{1});
                bonding[1] = new IntArrayList(new int[]{0,3, 2});
                bonding[2] = new IntArrayList(new int[]{1, 4});
                bonding[3] = new IntArrayList(new int[]{1});
                bonding[4] = new IntArrayList(new int[]{2, 5});
                bonding[5] = new IntArrayList(new int[]{4});



                double x3 = bondLengthCC - bondLengthCC * Math.cos(theta_CCC);
                double y3 = -bondLengthCC * Math.sin(theta_CCC);
                double x4 = bondLengthCC - bondLengthCC * Math.cos(theta_CCC);
                double y4 = bondLengthCC * Math.sin(theta_CCC);
                double xO = x3 + bondLengthCH3OH * Math.sin(theta_CCOH);
                double yO = y3 + bondLengthCH3OH * Math.cos(theta_CCOH);
                double xH = xO + bondLengthOH * Math.sin(theta_CH3OH);
                double yH = yO + bondLengthOH * Math.cos(theta_CH3OH);
                //Get Coordinates
                Vector3D posC1 = new Vector3D(new double[]{0, 0, 0});
                Vector3D posCH = new Vector3D(new double[]{bondLengthCC, 0, 0});
                Vector3D posCH2 = new Vector3D(new double[]{x3, y3, 0});
                Vector3D posC4 = new Vector3D(new double[]{x4, y4, 0});

                Vector3D posO = new Vector3D(new double[]{xO, yO, 0});
                Vector3D posH = new Vector3D(new double[]{xH, yH, 0});

                //CHECK types
                //Set Geometry
                species = new SpeciesBuilder(space)
                        .addAtom(typeCH3, posC1, "CH3")
                        .addAtom(typeCH, posCH, "CH")
                        .addAtom(typeCH2, posCH2, "CH3_3")
                        .addAtom(typeCH3, posC4, "CH3_")

                        .addAtom(typeO, posO, "O")
                        .addAtom(typeH, posH, "H")
                        .build();
            }



            else if (chemForm == ChemForm.benzene) {
                //TraPPE-UA
                //planar
                AtomType typeCH = new AtomType(Carbon.INSTANCE);

                atomTypes = new AtomType[]{typeCH};
                isFlex = false;
                //TraPPE Parameters
                double sigmaCH = 3.695; // Angstrom
                double epsilonCH = Kelvin.UNIT.toSim(50.5);
                double qCH = Electron.UNIT.toSim(0.0);
                double r = 1.40; //Angstrom

                //Construct Arrays
                sigma = new double[] {sigmaCH};
                epsilon = new double[] {epsilonCH};
                charge = new double[]{qCH};
                double x1 = Math.sqrt(3) * r / 2;
                double y1 = r/2;

                //Get Coordinates
                Vector3D posC1 = new Vector3D(new double[]{x1, y1, 0});
                Vector3D posC2 = new Vector3D(new double[]{x1, -y1, 0});
                Vector3D posC3 = new Vector3D(new double[]{0, -r, 0});
                Vector3D posC4 = new Vector3D(new double[]{-x1, -y1, 0});
                Vector3D posC5 = new Vector3D(new double[]{-x1, y1, 0});
                Vector3D posC6 = new Vector3D(new double[]{0, r, 0});


                //Set Geometry
                species = new SpeciesBuilder(space)
                        .addAtom(typeCH, posC1, "C1")
                        .addAtom(typeCH, posC2, "C2")
                        .addAtom(typeCH, posC3, "C3")
                        .addAtom(typeCH, posC4, "C4")
                        .addAtom(typeCH, posC5, "C5")
                        .addAtom(typeCH, posC6, "C6")

                        .build();
            }
            else if (chemForm == ChemForm.toluene) {
                //TraPPE 7 site UA
                //Atom in Compound
                AtomType typeC = new AtomType(Carbon.INSTANCE);
                AtomType typeCH = new AtomType(Carbon.INSTANCE);
                AtomType typeCH3 = new AtomType(Carbon.INSTANCE);

                atomTypes = new AtomType[]{typeC, typeCH, typeCH3};
                isFlex = false;
                //TraPPE Parameters
                double r = 1.40; // Angstrom
                double branch = 1.54; //Angstrom
                double theta = Degree.UNIT.toSim(30);
                double sigmaC = 3.88; // Angstrom
                double epsilonC = Kelvin.UNIT.toSim(21.0);
                double qC = Electron.UNIT.toSim(0);
                double sigmaCH = 3.695; // Angstrom
                double epsilonCH = Kelvin.UNIT.toSim(50.5);
                double qCH = Electron.UNIT.toSim(0);
                double sigmaCH3 = 3.75; // Angstrom
                double epsilonCH3 = Kelvin.UNIT.toSim(98);
                double qCH3 = Electron.UNIT.toSim(0);


                //Construct Arrays
                sigma = new double[] {sigmaC,sigmaCH, sigmaCH3};
                epsilon = new double[] {epsilonC, epsilonCH, epsilonCH3};
                charge = new double[]{qC, qCH, qCH3};

                double x1 = Math.sqrt(3) * r / 2;
                double y1 = r/2;

                        //Get Coordinates
                Vector3D posC1 = new Vector3D(new double[]{x1, y1, 0});
                Vector3D posC2 = new Vector3D(new double[]{x1, -y1, 0});
                Vector3D posC3 = new Vector3D(new double[]{0, -r, 0});
                Vector3D posC4 = new Vector3D(new double[]{-x1, -y1, 0});
                Vector3D posC5 = new Vector3D(new double[]{-x1, y1, 0});
                Vector3D posC6 = new Vector3D(new double[]{0, r, 0});
                Vector3D posC7 = new Vector3D(new double[]{x1 + branch * Math.cos(theta), y1 + branch * Math.sin(theta), 0});

                //Set Geometry
                species = new SpeciesBuilder(space)
                        .addAtom(typeC, posC1, "C1")
                        .addAtom(typeCH, posC2, "C2")
                        .addAtom(typeCH, posC3, "C3")
                        .addAtom(typeCH, posC4, "C4")
                        .addAtom(typeCH, posC5, "C5")
                        .addAtom(typeCH, posC6, "C6")
                        .addAtom(typeCH3, posC7, "C7")

                        .build();
            }
            else if (chemForm == ChemForm.oxylene) {
                //TraPPE UA
                //Atom in Compound
                AtomType typeC = new AtomType(Carbon.INSTANCE);
                AtomType typeCH = new AtomType(Carbon.INSTANCE);
                AtomType typeCH3 = new AtomType(Carbon.INSTANCE);

                atomTypes = new AtomType[]{typeC, typeCH, typeCH3};
                isFlex = false;
                //TraPPE Parameters
                double r = 1.40; // Angstrom
                double b = 1.54; //Angstrom
                double theta = Degree.UNIT.toSim(30);

                double sigmaC = 3.88; // Angstrom
                double epsilonC = Kelvin.UNIT.toSim(21.0);
                double qC = Electron.UNIT.toSim(0);
                double sigmaCH = 3.695; // Angstrom
                double epsilonCH = Kelvin.UNIT.toSim(50.5);
                double qCH = Electron.UNIT.toSim(0);
                double sigmaCH3 = 3.75; // Angstrom
                double epsilonCH3 = Kelvin.UNIT.toSim(98);
                double qCH3 = Electron.UNIT.toSim(0);


                //Construct Arrays
                sigma = new double[] {sigmaC,sigmaCH, sigmaCH3};
                epsilon = new double[] {epsilonC, epsilonCH, epsilonCH3};
                charge = new double[]{qC, qCH, qCH3};

                double x1 = Math.sqrt(3) * r / 2;
                double y1 = r/2;

                //Get Coordinates
                Vector3D posC1 = new Vector3D(new double[]{x1, y1, 0});
                Vector3D posC2 = new Vector3D(new double[]{x1, -y1, 0});
                Vector3D posC3 = new Vector3D(new double[]{0, -r, 0});
                Vector3D posC4 = new Vector3D(new double[]{-x1, -y1, 0});
                Vector3D posC5 = new Vector3D(new double[]{-x1, y1, 0});
                Vector3D posC6 = new Vector3D(new double[]{0, r, 0});
                Vector3D posC7 = new Vector3D(new double[]{x1 + b * Math.cos(theta), y1 + b * Math.sin(theta), 0});
                Vector3D posC8 = new Vector3D(new double[]{x1 + b * Math.cos(theta), -y1 - b * Math.sin(theta), 0});

                //Set Geometry
                species = new SpeciesBuilder(space)
                        .addAtom(typeC, posC1, "C1")
                        .addAtom(typeC, posC2, "C2")
                        .addAtom(typeCH, posC3, "C3")
                        .addAtom(typeCH, posC4, "C4")
                        .addAtom(typeCH, posC5, "C5")
                        .addAtom(typeCH, posC6, "C6")
                        .addAtom(typeCH3, posC7, "C7")
                        .addAtom(typeCH3, posC8, "C8")

                        .build();
            }
            else if (chemForm == ChemForm.pxylene) {
                //TraPPE UA
                //Atom in Compound
                AtomType typeC = new AtomType(Carbon.INSTANCE);
                AtomType typeCH = new AtomType(Carbon.INSTANCE);
                AtomType typeCH3 = new AtomType(Carbon.INSTANCE);

                atomTypes = new AtomType[]{typeC, typeCH, typeCH3};
                isFlex = false;
                //TraPPE Parameters
                double r = 1.40; // Angstrom
                double b = 1.54; //Angstrom
                double sigmaC = 3.88; // Angstrom
                double theta = Degree.UNIT.toSim(30);

                double epsilonC = Kelvin.UNIT.toSim(21.0);
                double qC = Electron.UNIT.toSim(0);
                double sigmaCH = 3.695; // Angstrom
                double epsilonCH = Kelvin.UNIT.toSim(50.5);
                double qCH = Electron.UNIT.toSim(0);
                double sigmaCH3 = 3.75; // Angstrom
                double epsilonCH3 = Kelvin.UNIT.toSim(98);
                double qCH3 = Electron.UNIT.toSim(0);


                //Construct Arrays
                sigma = new double[] {sigmaC,sigmaCH, sigmaCH3};
                epsilon = new double[] {epsilonC, epsilonCH, epsilonCH3};
                charge = new double[]{qC, qCH, qCH3};

                double x1 = Math.sqrt(3) * r / 2;
                double y1 = r/2;

                //Get Coordinates
                Vector3D posC1 = new Vector3D(new double[]{x1, y1, 0});
                Vector3D posC2 = new Vector3D(new double[]{x1, -y1, 0});
                Vector3D posC3 = new Vector3D(new double[]{0, -r, 0});
                Vector3D posC4 = new Vector3D(new double[]{-x1, -y1, 0});
                Vector3D posC5 = new Vector3D(new double[]{-x1, y1, 0});
                Vector3D posC6 = new Vector3D(new double[]{0, r, 0});
                Vector3D posC7 = new Vector3D(new double[]{x1 + b * Math.cos(theta), y1 + b * Math.sin(theta), 0});
                Vector3D posC8 = new Vector3D(new double[]{-x1 - b * Math.cos(theta), -y1 - b * Math.sin(theta), 0});
                //Set Geometry
                species = new SpeciesBuilder(space)
                        .addAtom(typeC, posC1, "C1")
                        .addAtom(typeC, posC4, "C4")
                        .addAtom(typeCH, posC3, "C3")
                        .addAtom(typeCH, posC2, "C2")
                        .addAtom(typeCH, posC5, "C5")
                        .addAtom(typeCH, posC6, "C6")
                        .addAtom(typeCH3, posC7, "C7")
                        .addAtom(typeCH3, posC8, "C8")

                        .build();
            }
            else if (chemForm == ChemForm.mxylene) {
                //TraPPE UA
                //Atom in Compound
                AtomType typeC = new AtomType(Carbon.INSTANCE);
                AtomType typeCH = new AtomType(Carbon.INSTANCE);
                AtomType typeCH3 = new AtomType(Carbon.INSTANCE);

                atomTypes = new AtomType[]{typeC, typeCH, typeCH3};
                isFlex = false;
                //TraPPE Parameters
                double r = 1.40; // Angstrom
                double b = 1.54; //Angstrom
                double sigmaC = 3.88; // Angstrom
                double theta = Degree.UNIT.toSim(30);

                double epsilonC = Kelvin.UNIT.toSim(21.0);
                double qC = Electron.UNIT.toSim(0);
                double sigmaCH = 3.695; // Angstrom
                double epsilonCH = Kelvin.UNIT.toSim(50.5);
                double qCH = Electron.UNIT.toSim(0);
                double sigmaCH3 = 3.75; // Angstrom
                double epsilonCH3 = Kelvin.UNIT.toSim(98);
                double qCH3 = Electron.UNIT.toSim(0);


                //Construct Arrays
                sigma = new double[] {sigmaC,sigmaCH, sigmaCH3};
                epsilon = new double[] {epsilonC, epsilonCH, epsilonCH3};
                charge = new double[]{qC, qCH, qCH3};

                double x1 = Math.sqrt(3) * r / 2;
                double y1 = r/2;

                //Get Coordinates
                Vector3D posC1 = new Vector3D(new double[]{x1, y1, 0});
                Vector3D posC2 = new Vector3D(new double[]{x1, -y1, 0});
                Vector3D posC3 = new Vector3D(new double[]{0, -r, 0});
                Vector3D posC4 = new Vector3D(new double[]{-x1, -y1, 0});
                Vector3D posC5 = new Vector3D(new double[]{-x1, y1, 0});
                Vector3D posC6 = new Vector3D(new double[]{0, r, 0});
                Vector3D posC7 = new Vector3D(new double[]{x1 + b * Math.cos(theta), y1 + b * Math.sin(theta), 0});
                Vector3D posC8 = new Vector3D(new double[]{0, -r - b, 0});

                //Set Geometry
                species = new SpeciesBuilder(space)
                        .addAtom(typeC, posC1, "C1")
                        .addAtom(typeC, posC3, "C3")
                        .addAtom(typeCH, posC4, "C4")
                        .addAtom(typeCH, posC2, "C2")
                        .addAtom(typeCH, posC5, "C5")
                        .addAtom(typeCH, posC6, "C6")
                        .addAtom(typeCH3, posC7, "C7")
                        .addAtom(typeCH3, posC8, "C8")

                        .build();
            }

            else if (chemForm == ChemForm.ethylbenzene) {
                //TraPPE UA
                // ok
                AtomType typeC = new AtomType(Carbon.INSTANCE);
                AtomType typeCH = new AtomType(Carbon.INSTANCE);
                AtomType typeCH2 = new AtomType(Carbon.INSTANCE);
                AtomType typeCH3 = new AtomType(Carbon.INSTANCE);

                atomTypes = new AtomType[]{typeC, typeCH, typeCH2, typeCH3};
                isFlex = true;
                //TraPPE Parameters
                double r = 1.40; // Angstrom
                double b = 1.54; //Angstrom
                double thetaCCCy = Degree.UNIT.toSim(114);
                double theta = Degree.UNIT.toSim(36);

                double sigmaC = 3.88; // Angstrom
                double epsilonC = Kelvin.UNIT.toSim(21.0);
                double qC = Electron.UNIT.toSim(0);
                double sigmaCH = 3.695; // Angstrom
                double epsilonCH = Kelvin.UNIT.toSim(50.5);
                double qCH = Electron.UNIT.toSim(0);
                double sigmaCH2 = 3.95; // Angstrom
                double epsilonCH2 = Kelvin.UNIT.toSim(46);
                double qCH2 = Electron.UNIT.toSim(0);
                double sigmaCH3 = 3.75; // Angstrom
                double epsilonCH3 = Kelvin.UNIT.toSim(98);
                double qCH3 = Electron.UNIT.toSim(0);
                double k_thetaCCCy = Kelvin.UNIT.toSim(62500);


                //Construct Arrays
                sigma = new double[] {sigmaC,sigmaCH, sigmaCH2, sigmaCH3};
                epsilon = new double[] {epsilonC, epsilonCH, epsilonCH2, epsilonCH3};
                charge = new double[]{qC, qCH, qCH2, qCH3};
                theta_eq = new double[]{thetaCCCy};
                k_theta = new double[]{k_thetaCCCy};

                double x1 = Math.sqrt(3) * r / 2;
                double y1 = r/2;
                double x7 = x1 + b * Math.cos(Math.PI/6);
                double y7 = y1 + b * Math.sin(Math.PI/6);
                double e0 = Kelvin.UNIT.toSim(131);
                double e1 = Degree.UNIT.toSim(180);
                a = new double[][]{{e0, e1}};

                triplets = new int[][]{{0, 6, 7}};
                quads = new int[][]{{1, 0, 6, 7}, {5, 0, 6, 7}};
                bonding = new IntArrayList[9];
                bonding[0] = new IntArrayList(new int[]{1, 5, 6});
                bonding[1] = new IntArrayList(new int[]{0,2});
                bonding[2] = new IntArrayList(new int[]{1, 3});
                bonding[3] = new IntArrayList(new int[]{2, 4});
                bonding[4] = new IntArrayList(new int[]{3, 5});
                bonding[5] = new IntArrayList(new int[]{0, 4});
                bonding[6] = new IntArrayList(new int[]{0, 7});
                bonding[7] = new IntArrayList(new int[]{6});

                //Get Coordinates
                Vector3D posC1 = new Vector3D(new double[]{x1, y1, 0});
                Vector3D posC2 = new Vector3D(new double[]{x1, -y1, 0});
                Vector3D posC3 = new Vector3D(new double[]{0, -r, 0});
                Vector3D posC4 = new Vector3D(new double[]{-x1, -y1, 0});
                Vector3D posC5 = new Vector3D(new double[]{-x1, y1, 0});
                Vector3D posC6 = new Vector3D(new double[]{0, r, 0});
                Vector3D posC7 = new Vector3D(new double[]{x7, y7, 0});
                Vector3D posC8 = new Vector3D(new double[]{x7 + b * Math.cos(theta), y7 - b * Math.sin(theta), 0});


                //Set Geometry
                species = new SpeciesBuilder(space)
                        .addAtom(typeC, posC1, "C1")
                        .addAtom(typeCH, posC2, "C2")
                        .addAtom(typeCH, posC3, "C3")
                        .addAtom(typeCH, posC4, "C4")
                        .addAtom(typeCH, posC5, "C5")
                        .addAtom(typeCH, posC6, "C6")
                        .addAtom(typeCH2, posC7, "C7")
                        .addAtom(typeCH3, posC8, "C8")

                        .build();
            }

//            else if (chemForm == ChemForm.ethane) {
//                //TraPPE-EH
//                //Atom in Compound
//                //avogadro
//                AtomType typeC = new AtomType(Carbon.INSTANCE);
//                AtomType typeM = new AtomType(elementM);
//
//                atomTypes = new AtomType[]{typeC,typeM};
//                isFlex = false;
//                //TraPPE Parameters
//                double bondLengthCM = 0.55; // Angstrom
//                double bondLengthCC = 1.5350; //Angstrom
//                double thetaCCH = Degree.UNIT.toSim(110.70);
//                double a120 = Degree.UNIT.toSim(120);
//                double a240 = Degree.UNIT.toSim(240);
//                double sigmaC = 3.30; // Angstrom
//                double epsilonC = Kelvin.UNIT.toSim(4.00);
//                double qC = Electron.UNIT.toSim(0.0);
//                double sigmaM = 3.31; // Angstrom
//                double epsilonM = Kelvin.UNIT.toSim(15.30);
//                double qM = Electron.UNIT.toSim(0.000);
//                double a01 = Kelvin.UNIT.toSim(717);
//                //Construct Arrays
//                sigma = new double[] {sigmaC,sigmaM};
//                epsilon = new double[] {epsilonC,epsilonM};
//                charge = new double[]{qC, qM};
//                a = new double[][]{{0, a01, 0, 0}};
//
//
//                //Get Coordinates
//                Vector3D posC1 = new Vector3D(new double[]{0.7674410120,     -0.0579174161  ,   -0.0637971165});
//                Vector3D posC2 = new Vector3D(new double[]{-0.7632720540,      0.0567215202,     -0.0644482175});
//                Vector3D posM1 = new Vector3D(new double[]{0.9474423656,     -0.3387486927,      0.3735051894});
//                Vector3D posM2 = new Vector3D(new double[]{0.9394900365 ,    -0.3210710943,     -0.5150722819});
//                Vector3D posM3 = new Vector3D(new double[]{0.9997317223,      0.4406138058,     -0.0610318522});
//                Vector3D posM4 = new Vector3D(new double[]{-0.9436953068,      0.3368882196,      0.3731064310});
//                Vector3D posM5 = new Vector3D(new double[]{-0.9348859373 ,     0.3205606343,     -0.5154887887});
//                Vector3D posM6 = new Vector3D(new double[]{-0.9955649581,     -0.4418131574,     -0.0626644104});
//
//                //Set Geometry
//                species = new SpeciesBuilder(space)
//                        .addAtom(typeC, posC1, "C1")
//                        .addAtom(typeC, posC2, "C2")
//                        .addAtom(typeM, posM1, "M1")
//                        .addAtom(typeM, posM2, "M2")
//                        .addAtom(typeM, posM3, "M3")
//                        .addAtom(typeM, posM4, "M4")
//                        .addAtom(typeM, posM5, "M5")
//                        .addAtom(typeM, posM6, "M6")
//
//                        .build();
//            }
//            else if (chemForm == ChemForm.propane) {
//                //TraPPE-EH
//                //Atom in Compound
//                //Avogadro
//                AtomType typeC = new AtomType(Carbon.INSTANCE);
//                AtomType typemidC = new AtomType(Carbon.INSTANCE);
//                AtomType typeM = new AtomType(elementM);
//
//                atomTypes = new AtomType[]{typeC, typemidC, typeM};
//                isFlex = true;
//                //TraPPE Parameters
//                double bondLengthCM = 0.55; // Angstrom
//                double bondLengthCC = 1.5350; //Angstrom
//                double alpha = Degree.UNIT.toSim(53.90);
//                double thetaCCH = Degree.UNIT.toSim(110.70);
//                double thetaCCC = Degree.UNIT.toSim(112.70);
//                double a120 = Degree.UNIT.toSim(120);
//                double a240 = Degree.UNIT.toSim(240);
//                double sigmaC = 3.30; // Angstrom
//                double epsilonC = Kelvin.UNIT.toSim(4.00);
//                double qC = Electron.UNIT.toSim(0.0);
//                double sigmamidC = 3.65; // Angstrom
//                double epsilonmidC = Kelvin.UNIT.toSim(5.00);
//                double qmidC = Electron.UNIT.toSim(0.0);
//                double sigmaM = 3.31; // Angstrom
//                double epsilonM = Kelvin.UNIT.toSim(15.30);
//                double qM = Electron.UNIT.toSim(0.000);
//                double k = Kelvin.UNIT.toSim(58765);
//                double a01 = Kelvin.UNIT.toSim(854);
//                //Construct Arrays
//                sigma = new double[] {sigmaC,sigmamidC, sigmaM};
//                epsilon = new double[] {epsilonC,epsilonmidC, epsilonM};
//                charge = new double[]{qC, qmidC, qM};
//                k_theta = new double[]{k};
//                theta_eq = new double[]{thetaCCC};
//                a = new double[][]{{0, a01, 0, 0}}; //fix, 3-1-2-4
//
//
//                double x3 = bondLengthCC - bondLengthCC * Math.cos(thetaCCC);
//                double y3 = bondLengthCC * Math.sin(thetaCCC);
//                //Get Coordinates
//                Vector3D posC1 = new Vector3D(new double[]{1.2520730479,     -0.2644906931,     -0.1533486375});
//                Vector3D posC2 = new Vector3D(new double[]{0.0099234689,      0.6324679426,     -0.0597824612});
//                Vector3D posC3 = new Vector3D(new double[]{-1.2958658885 ,    -0.1704726300,      0.0203773238});
//                Vector3D posM1 = new Vector3D(new double[]{1.2785732026,     -0.6058207170,      0.2771064201});
//                Vector3D posM2 = new Vector3D(new double[]{1.2323681877,     -0.5707199860,     -0.6097862988});
//                Vector3D posM3 = new Vector3D(new double[]{1.7030456657 ,     0.0491621611,     -0.1806540407});
//                Vector3D posM4 = new Vector3D(new double[]{0.0613163156,      0.9576922390,      0.3807716584});
//                Vector3D posM5 = new Vector3D(new double[]{0.0003552066,      0.9539846380,     -0.5059163682});
//                Vector3D posM6 = new Vector3D(new double[]{-1.3069583182,     -0.4786267642,      0.4758089735});
//                Vector3D posM7 = new Vector3D(new double[]{-1.7235437146 ,     0.1747612163,      0.0405086015});
//                Vector3D posM8 = new Vector3D(new double[]{-1.3538001794,     -0.4839632176,     -0.4278048595});
//
//                //Set Geometry
//                species = new SpeciesBuilder(space)
//                        .addAtom(typeC, posC1, "C1")
//                        .addAtom(typemidC, posC2, "C2")
//                        .addAtom(typeC, posC3, "C3")
//                        .addAtom(typeM, posM1, "M1")
//                        .addAtom(typeM, posM2, "M2")
//                        .addAtom(typeM, posM3, "M3")
//                        .addAtom(typeM, posM4, "M4")
//                        .addAtom(typeM, posM5, "M5")
//                        .addAtom(typeM, posM6, "M6")
//                        .addAtom(typeM, posM7, "M7")
//                        .addAtom(typeM, posM8, "M8")
//                        .build();
//                System.out.println("Carbon Positions: " + Arrays.toString(new Vector3D[]{posC1, posC2, posC3}));
//                System.out.println("Hydrogen Positions: " + Arrays.toString(new Vector3D[]{posM1, posM2, posM3, posM4, posM5, posM6, posM7, posM8}));
//            }
            else if (chemForm == ChemForm.ethane) {
                //TraPPE-UA
                //Atom in Compound
                //QC checked
                AtomType typeCH3 = new AtomType(Carbon.INSTANCE);

                atomTypes = new AtomType[]{typeCH3};
                isFlex = false;
                //TraPPE Parameters
                double bondLengthCHxCHy = 1.54; // Angstrom

                double sigmaC = 3.75; // Angstrom
                double epsilonC = Kelvin.UNIT.toSim(98);
                double qC = Electron.UNIT.toSim(0.0);

                //Construct Arrays
                sigma = new double[] {sigmaC};
                epsilon = new double[] {epsilonC};
                charge = new double[]{qC};

                //Get Coordinates
                Vector3D posC1 = new Vector3D(new double[]{0, 0, 0});
                Vector3D posC2 = new Vector3D(new double[]{bondLengthCHxCHy, 0, 0});

                //Set Geometry
                species = new SpeciesBuilder(space)
                        .addAtom(typeCH3, posC1, "C1")
                        .addAtom(typeCH3, posC2, "C2")
                        .build();
            }
            else if (chemForm == ChemForm.propane) {
                //TraPPE-UA
                //Atom in Compound
                //Avogadro
                AtomType typeCH3 = new AtomType(Carbon.INSTANCE);
                AtomType typeCH2 = new AtomType(Carbon.INSTANCE);

                atomTypes = new AtomType[]{typeCH3, typeCH2};
                isFlex = true;
                //TraPPE Parameters
                double bondLengthCHxCHy = 1.54; // Angstrom
                double thetaCCH = Degree.UNIT.toSim(114);
                double sigmaCH3 = 3.75; // Angstrom
                double epsilonCH3 = Kelvin.UNIT.toSim(98);
                double qCH3 = Electron.UNIT.toSim(0.0);
                double sigmaCH2 = 3.95; // Angstrom
                double epsilonCH2 = Kelvin.UNIT.toSim(46);
                double qCH2 = Electron.UNIT.toSim(0.0);
                double kCCC = Kelvin.UNIT.toSim(62500);
                double yy = bondLengthCHxCHy*Math.sin(thetaCCH) ;
                //Construct Arrays
                sigma = new double[] {sigmaCH3, sigmaCH2};
                epsilon = new double[] {epsilonCH3, epsilonCH2};
                charge = new double[]{qCH3, qCH2};
                k_theta = new double[]{kCCC};
                theta_eq = new double[]{thetaCCH};
                triplets = new int[][]{{0, 1, 2}};
                bonding = new IntArrayList[3];
                bonding[0] = new IntArrayList(new int[]{1});
                bonding[1] = new IntArrayList(new int[]{0,2});
                bonding[2] = new IntArrayList(new int[]{1});


                //Get Coordinates
                Vector3D posC1 = new Vector3D(new double[]{0, -yy/3, 0});
                Vector3D posC2 = new Vector3D(new double[]{bondLengthCHxCHy, -yy/3, 0});
                Vector3D posC3 = new Vector3D(new double[]{bondLengthCHxCHy - bondLengthCHxCHy * Math.cos(thetaCCH), 2*yy/3, 0});
                double posCavg = (2 * bondLengthCHxCHy  - bondLengthCHxCHy * Math.cos(thetaCCH))/3;
                posC1 = new Vector3D(new double[]{0 - posCavg, -yy / 3, 0});
                posC2 = new Vector3D(new double[]{bondLengthCHxCHy - posCavg, -yy / 3, 0});
                posC3 = new Vector3D(new double[]{bondLengthCHxCHy - bondLengthCHxCHy * Math.cos(thetaCCH) - posCavg, 2 * yy / 3, 0});

                //Set Geometry
                species = new SpeciesBuilder(space)
                        .addAtom(typeCH3, posC1, "C1")
                        .addAtom(typeCH2, posC2, "C2")
                        .addAtom(typeCH3, posC3, "C3")
                        .build();

            }
            else if (chemForm == ChemForm.butane) {
                //TraPPE-UA
                //Atom in Compound
                //avogadro trans
                AtomType typeCH3 = new AtomType(Carbon.INSTANCE);
                AtomType typeCH2 = new AtomType(Carbon.INSTANCE);

                atomTypes = new AtomType[]{typeCH3, typeCH2};
                isFlex = true;
                //TraPPE Parameters
                double bondLengthCHxCHy = 1.54; // Angstrom
                double thetaCCH = Degree.UNIT.toSim(114);

                double sigmaCH3 = 3.75; // Angstrom
                double epsilonCH3 = Kelvin.UNIT.toSim(98);
                double qCH3 = Electron.UNIT.toSim(0.0);
                double sigmaCH2 = 3.95; // Angstrom
                double epsilonCH2 = Kelvin.UNIT.toSim(46);
                double qCH2 = Electron.UNIT.toSim(0.0);
                double k = Kelvin.UNIT.toSim(62500);
                double a00 = 0;
                double a01 = 355.03;
                double a02 = -68.19;
                double a03 = 791.32;
                //Construct Arrays
                sigma = new double[] {sigmaCH3, sigmaCH2};
                epsilon = new double[] {epsilonCH3, epsilonCH2};
                charge = new double[]{qCH3, qCH2};
                k_theta = new double[]{k, k};
                theta_eq = new double[]{thetaCCH, thetaCCH};
                a = new double[][]{{a00, a01, a02, a03}};
                triplets = new int[][]{{0, 1, 2}, {1, 2, 3}};
                quads = new int[][]{{0, 1, 2, 3}};
                bonding = new IntArrayList[4];
                bonding[0] = new IntArrayList(new int[]{1});
                bonding[1] = new IntArrayList(new int[]{0,2});
                bonding[2] = new IntArrayList(new int[]{1, 3});
                bonding[3] = new IntArrayList(new int[]{2});


                //Get Coordinates

                double m3 = typeCH3.getMass(), m2 = typeCH2.getMass();
                double x1 = -bondLengthCHxCHy, x2 = 0, x3 = -bondLengthCHxCHy*Math.cos(thetaCCH);
                double x4 = x3 + bondLengthCHxCHy;
                double y1 = 0, y2 = 0, y3 = bondLengthCHxCHy*Math.sin(thetaCCH);
                double y4 = y3;
                double mm = 2 * m3 + 2 * m2;
                double xCOM = (m3*x1 + m2*x2 + m2*x3 + m3*x4) / mm;
                double yCOM = (m3*y1 + m2*y2 + m2*y3 + m3*y4) / mm;

                Vector3D posC1 = new Vector3D(new double[]{x1-xCOM, y1-yCOM, 0});
                Vector3D posC2 = new Vector3D(new double[]{x2-xCOM, y2-yCOM, 0});
                Vector3D posC3 = new Vector3D(new double[]{x3-xCOM, y3-yCOM, 0});
                Vector3D posC4 = new Vector3D(new double[]{x4-xCOM, y4-yCOM, 0});

                System.out.println("Carbon Positions: " + Arrays.toString(new Vector3D[]{posC1, posC2, posC3, posC4}));

                //Set Geometry
                species = new SpeciesBuilder(space)
                        .addAtom(typeCH3, posC1, "C1")
                        .addAtom(typeCH2, posC2, "C2")
                        .addAtom(typeCH2, posC3, "C3")
                        .addAtom(typeCH3, posC4, "C4")

                        .build();
            }
            else if (chemForm == ChemForm.methane) {
                //TraPPE-UA
                //Atom in Compound
                //avogadro not needed
                AtomType typeCH4 = new AtomType(Carbon.INSTANCE);

                atomTypes = new AtomType[]{typeCH4};
                isFlex = false;
                //TraPPE Parameters
                double sigmaCH4 = 3.73; // Angstrom
                double epsilonCH4 = Kelvin.UNIT.toSim(148);
                double qCH4 = Electron.UNIT.toSim(0.0);

                //Construct Arrays
                sigma = new double[] {sigmaCH4};
                epsilon = new double[] {epsilonCH4};
                charge = new double[]{qCH4};

                //Get Coordinates
                Vector3D posC1 = new Vector3D(new double[]{0, 0, 0});

                //Set Geometry
                species = new SpeciesBuilder(space)
                        .addAtom(typeCH4, posC1, "C1")

                        .build();
            }
            else if (chemForm == ChemForm.ethene) {
                //TraPPE-UA
                //Atom in Compound
                //avogadro not needed
                AtomType typeCH2 = new AtomType(Carbon.INSTANCE);
                isFlex = false;
                atomTypes = new AtomType[]{typeCH2};

                //TraPPE Parameters
                double bondLengthCC = 1.33; //Angstrom

                double sigmaCH2 = 3.675; // Angstrom
                double epsilonCH2 = Kelvin.UNIT.toSim(85.00);
                double qCH2 = Electron.UNIT.toSim(0.0);

                //Construct Arrays
                sigma = new double[] {sigmaCH2};
                epsilon = new double[] {epsilonCH2};
                charge = new double[]{qCH2};

                //Get Coordinates
                Vector3D posC1 = new Vector3D(new double[]{0, 0, 0});
                Vector3D posC2 = new Vector3D(new double[]{bondLengthCC, 0, 0});

                //Set Geometry
                species = new SpeciesBuilder(space)
                        .addAtom(typeCH2, posC1, "C1")
                        .addAtom(typeCH2, posC2, "C2")

                        .build();
            }
            else if (chemForm == ChemForm.propene) {
                //TraPPE-UA
                //Atom in Compound
                //checked
                AtomType typeCH3 = new AtomType(Carbon.INSTANCE);
                AtomType typeCH = new AtomType(Carbon.INSTANCE);
                AtomType typeCH2 = new AtomType(Carbon.INSTANCE);

                atomTypes = new AtomType[]{typeCH3, typeCH, typeCH2};
                isFlex = true;
                //TraPPE Parameters
                double bondLengthCHxCHy = 1.54; // Angstrom
                double bondLengthCHxChy_double = 1.33; //Angstrom
                double thetaeq = Degree.UNIT.toSim(119.70);
                double sigmaCH3 = 3.75; // Angstrom
                double epsilonCH3 = Kelvin.UNIT.toSim(98.00);
                double qCH3 = Electron.UNIT.toSim(0.0);
                double sigmaCH = 3.73; // Angstrom
                double epsilonCH = Kelvin.UNIT.toSim(47.00);
                double qCH = Electron.UNIT.toSim(0.0);
                double sigmaCH2 = 3.675; // Angstrom
                double epsilonCH2 = Kelvin.UNIT.toSim(85.00);
                double qCH2 = Electron.UNIT.toSim(0.000);
                double k = Kelvin.UNIT.toSim(70420);

                //Construct Arrays
                sigma = new double[] {sigmaCH3,sigmaCH, sigmaCH2};
                epsilon = new double[] {epsilonCH3,epsilonCH, epsilonCH2};
                charge = new double[]{qCH3, qCH, qCH2};
                theta_eq = new double[]{thetaeq};
                k_theta = new double[]{k};
                triplets = new int[][]{{0, 1, 2}};
                bonding = new IntArrayList[3];
                bonding[0] = new IntArrayList(new int[]{1});
                bonding[1] = new IntArrayList(new int[]{0,2});
                bonding[2] = new IntArrayList(new int[]{1});

                //Get Coordinates

                Vector3D posCH3 = new Vector3D(new double[]{0, 0, 0});
                Vector3D posCH = new Vector3D(new double[]{bondLengthCHxCHy, 0, 0});
                Vector3D posCH2 = new Vector3D(new double[]{bondLengthCHxCHy - bondLengthCHxCHy*Math.cos(thetaeq), bondLengthCHxCHy*Math.sin(thetaeq), 0});

                System.out.println(posCH2);

                //Set Geometry
                species = new SpeciesBuilder(space)
                        .addAtom(typeCH3, posCH3, "CH3")
                        .addAtom(typeCH, posCH, "CH")
                        .addAtom(typeCH2, posCH2, "CH2")
                        .build();
            }
            else if (chemForm == ChemForm.butadiene) {
                //TraPPE-UA
                //Atom in Compound
                //checked
                AtomType typeCH = new AtomType(Carbon.INSTANCE);
                AtomType typeCH2 = new AtomType(Carbon.INSTANCE);

                atomTypes = new AtomType[]{typeCH, typeCH2};
                isFlex = true;
                //TraPPE Parameters
                double bondLengthCHxCHy = 1.54; // Angstrom
                double bondLengthCHxChy_double = 1.33; //Angstrom

                double thetaeq = Degree.UNIT.toSim(119.70);
                double sigmaCH = 3.710; // Angstrom
                double epsilonCH = Kelvin.UNIT.toSim(52.0);
                double qCH = Electron.UNIT.toSim(0.0);
                double sigmaCH2 = 3.675; // Angstrom
                double epsilonCH2 = Kelvin.UNIT.toSim(85.00);
                double qCH2 = Electron.UNIT.toSim(0.000);
                double k = Kelvin.UNIT.toSim(70420);
                double a0prime = Kelvin.UNIT.toSim(2034.58);
                double a1prime = Kelvin.UNIT.toSim(531.57);
                double a2prime = Kelvin.UNIT.toSim(-1239.35);
                double a3prime = Kelvin.UNIT.toSim(460.04);
                double a00 = a0prime - a1prime + a2prime - a3prime;
                //Construct Arrays
                sigma = new double[] {sigmaCH, sigmaCH2};
                epsilon = new double[] {epsilonCH, epsilonCH2};
                charge = new double[]{qCH, qCH2};
                theta_eq = new double[]{thetaeq, thetaeq}; //123, 234
                k_theta = new double[]{k, k}; //123, 234
                a = new double[][]{{a00, a1prime, -a2prime, a3prime}};
                triplets = new int[][]{{0, 1, 2}, {1, 2, 3}};
                quads = new int[][]{{0, 1, 2, 3}};
                bonding = new IntArrayList[4];
                bonding[0] = new IntArrayList(new int[]{1});
                bonding[1] = new IntArrayList(new int[]{0,2});
                bonding[2] = new IntArrayList(new int[]{1, 3});
                bonding[3] = new IntArrayList(new int[]{2});

                double x3 = - bondLengthCHxCHy * Math.cos(thetaeq) + bondLengthCHxChy_double;
                double y3 = bondLengthCHxCHy * Math.sin(thetaeq);

                //Get Coordinates
                Vector3D posC1 = new Vector3D(new double[]{0, 0, 0});
                Vector3D posC2 = new Vector3D(new double[]{bondLengthCHxChy_double, 0, 0});
                Vector3D posC3 = new Vector3D(new double[]{x3, y3,   0});
                Vector3D posC4 = new Vector3D(new double[]{x3 + bondLengthCHxChy_double, y3, 0});

                System.out.println("posC1: " + posC1 + "posC2:" + posC2 + "posC3:" + posC3 + "posC4:" + posC4);

                //Set Geometry
                species = new SpeciesBuilder(space)
                        .addAtom(typeCH2, posC1, "C1")
                        .addAtom(typeCH, posC2, "C2")
                        .addAtom(typeCH, posC3, "C3")
                        .addAtom(typeCH2, posC4, "C4")

                        .build();
            }

            else {
                throw new RuntimeException("unrecognized chem form");
            }

        }



        public void buildPotentials(SpeciesManager sm) {
            if(chemForm == ChemForm.N2) {
                P2PotentialGroupBuilder.ModelParams modelParams = new P2PotentialGroupBuilder.ModelParams(atomTypes, sigma, epsilon, charge);
                potentialGroup = P2PotentialGroupBuilder.P2PotentialGroupBuilder(space, sm, modelParams, null);
            }
            else if (chemForm == ChemForm.O2) {
                P2PotentialGroupBuilder.ModelParams modelParams = new P2PotentialGroupBuilder.ModelParams(atomTypes, sigma, epsilon, charge);
                potentialGroup = P2PotentialGroupBuilder.P2PotentialGroupBuilder(space, sm, modelParams, null);
            }
            else if (chemForm == ChemForm.CO2) {
                P2PotentialGroupBuilder.ModelParams modelParams = new P2PotentialGroupBuilder.ModelParams(atomTypes, sigma, epsilon, charge);
                potentialGroup = P2PotentialGroupBuilder.P2PotentialGroupBuilder(space, sm, modelParams, null);
            }
            else if (chemForm == ChemForm.NH3) {
                P2PotentialGroupBuilder.ModelParams modelParams = new P2PotentialGroupBuilder.ModelParams(atomTypes, sigma, epsilon, charge);
                potentialGroup = P2PotentialGroupBuilder.P2PotentialGroupBuilder(space, sm, modelParams, null);
            }
            else if (chemForm == ChemForm.CH4) {
                P2PotentialGroupBuilder.ModelParams modelParams = new P2PotentialGroupBuilder.ModelParams(atomTypes, sigma, epsilon, charge);
                potentialGroup = P2PotentialGroupBuilder.P2PotentialGroupBuilder(space, sm, modelParams, null);
            }
            else if (chemForm == ChemForm.CH3OH) {
                P2PotentialGroupBuilder.ModelParams modelParams = new P2PotentialGroupBuilder.ModelParams(atomTypes, sigma, epsilon, charge);
                potentialGroup = P2PotentialGroupBuilder.P2PotentialGroupBuilder(space, sm, modelParams, null);
            }
            else if (chemForm == ChemForm.benzene) {
                P2PotentialGroupBuilder.ModelParams modelParams = new P2PotentialGroupBuilder.ModelParams(atomTypes, sigma, epsilon, charge);
                potentialGroup = P2PotentialGroupBuilder.P2PotentialGroupBuilder(space, sm, modelParams, null);
            }
            else if (chemForm == ChemForm.ethane) {
                P2PotentialGroupBuilder.ModelParams modelParams = new P2PotentialGroupBuilder.ModelParams(atomTypes, sigma, epsilon, charge);
                potentialGroup = P2PotentialGroupBuilder.P2PotentialGroupBuilder(space, sm, modelParams, null);
            }
            else if (chemForm == ChemForm.propane) {
                P2PotentialGroupBuilder.ModelParams modelParams = new P2PotentialGroupBuilder.ModelParams(atomTypes, sigma, epsilon, charge);
                potentialGroup = P2PotentialGroupBuilder.P2PotentialGroupBuilder(space, sm, modelParams, null);
            }
            else if (chemForm == ChemForm.butane) {
                P2PotentialGroupBuilder.ModelParams modelParams = new P2PotentialGroupBuilder.ModelParams(atomTypes, sigma, epsilon, charge);
                potentialGroup = P2PotentialGroupBuilder.P2PotentialGroupBuilder(space, sm, modelParams, null);
            }
            else if (chemForm == ChemForm.methane) {
                P2PotentialGroupBuilder.ModelParams modelParams = new P2PotentialGroupBuilder.ModelParams(atomTypes, sigma, epsilon, charge);
                potentialGroup = P2PotentialGroupBuilder.P2PotentialGroupBuilder(space, sm, modelParams, null);
            }
            else if (chemForm == ChemForm.ethene) {
                P2PotentialGroupBuilder.ModelParams modelParams = new P2PotentialGroupBuilder.ModelParams(atomTypes, sigma, epsilon, charge);
                potentialGroup = P2PotentialGroupBuilder.P2PotentialGroupBuilder(space, sm, modelParams, null);
            }
            else if (chemForm == ChemForm.propene) {
                P2PotentialGroupBuilder.ModelParams modelParams = new P2PotentialGroupBuilder.ModelParams(atomTypes, sigma, epsilon, charge);
                potentialGroup = P2PotentialGroupBuilder.P2PotentialGroupBuilder(space, sm, modelParams, null);
            }
            else if (chemForm == ChemForm.butadiene) {
                P2PotentialGroupBuilder.ModelParams modelParams = new P2PotentialGroupBuilder.ModelParams(atomTypes, sigma, epsilon, charge);
                potentialGroup = P2PotentialGroupBuilder.P2PotentialGroupBuilder(space, sm, modelParams, null);
            }
            else if (chemForm == ChemForm.ethanol) {
                P2PotentialGroupBuilder.ModelParams modelParams = new P2PotentialGroupBuilder.ModelParams(atomTypes, sigma, epsilon, charge);
                potentialGroup = P2PotentialGroupBuilder.P2PotentialGroupBuilder(space, sm, modelParams, null);
            }
            else if (chemForm == ChemForm.propan1ol) {
                P2PotentialGroupBuilder.ModelParams modelParams = new P2PotentialGroupBuilder.ModelParams(atomTypes, sigma, epsilon, charge);
                potentialGroup = P2PotentialGroupBuilder.P2PotentialGroupBuilder(space, sm, modelParams, null);
            }
            else if (chemForm == ChemForm.propan2ol) {
                P2PotentialGroupBuilder.ModelParams modelParams = new P2PotentialGroupBuilder.ModelParams(atomTypes, sigma, epsilon, charge);
                potentialGroup = P2PotentialGroupBuilder.P2PotentialGroupBuilder(space, sm, modelParams, null);
            }
            else if (chemForm == ChemForm.isobutanol) {
                P2PotentialGroupBuilder.ModelParams modelParams = new P2PotentialGroupBuilder.ModelParams(atomTypes, sigma, epsilon, charge);
                potentialGroup = P2PotentialGroupBuilder.P2PotentialGroupBuilder(space, sm, modelParams, null);
            }
            else if (chemForm == ChemForm.toluene) {
                P2PotentialGroupBuilder.ModelParams modelParams = new P2PotentialGroupBuilder.ModelParams(atomTypes, sigma, epsilon, charge);
                potentialGroup = P2PotentialGroupBuilder.P2PotentialGroupBuilder(space, sm, modelParams, null);
            }
            else if (chemForm == ChemForm.ethylbenzene) {
                P2PotentialGroupBuilder.ModelParams modelParams = new P2PotentialGroupBuilder.ModelParams(atomTypes, sigma, epsilon, charge);
                potentialGroup = P2PotentialGroupBuilder.P2PotentialGroupBuilder(space, sm, modelParams, null);
            }
            else if (chemForm == ChemForm.oxylene) {
                P2PotentialGroupBuilder.ModelParams modelParams = new P2PotentialGroupBuilder.ModelParams(atomTypes, sigma, epsilon, charge);
                potentialGroup = P2PotentialGroupBuilder.P2PotentialGroupBuilder(space, sm, modelParams, null);
            }
            else if (chemForm == ChemForm.mxylene) {
                P2PotentialGroupBuilder.ModelParams modelParams = new P2PotentialGroupBuilder.ModelParams(atomTypes, sigma, epsilon, charge);
                potentialGroup = P2PotentialGroupBuilder.P2PotentialGroupBuilder(space, sm, modelParams, null);
            }
            else if (chemForm == ChemForm.pxylene) {
                P2PotentialGroupBuilder.ModelParams modelParams = new P2PotentialGroupBuilder.ModelParams(atomTypes, sigma, epsilon, charge);
                potentialGroup = P2PotentialGroupBuilder.P2PotentialGroupBuilder(space, sm, modelParams, null);
            }
            else if (chemForm == ChemForm.water) {
                P2PotentialGroupBuilder.ModelParams modelParams = new P2PotentialGroupBuilder.ModelParams(atomTypes, sigma, epsilon, charge);
                potentialGroup = P2PotentialGroupBuilder.P2PotentialGroupBuilder(space, sm, modelParams, null);
            }



        }

    }

}
