/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.chem.elements.*;
import etomica.graph.model.Graph;
import etomica.graphics.ColorSchemeRandomByMolecule;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayBoxCanvasG3DSys;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.mcmove.MCMove;
import etomica.integrator.mcmove.MCMoveBox;
import etomica.math.SpecialFunctions;
import etomica.molecule.IMoleculeList;
import etomica.molecule.MoleculePositionCOM;
import etomica.potential.*;
import etomica.space.Space;
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
import etomica.virial.*;
import etomica.virial.cluster.*;
import etomica.virial.mcmove.*;
import etomica.virial.wheatley.ClusterWheatleyHS;
import etomica.virial.wheatley.ClusterWheatleySoftDerivatives;
import etomica.virial.wheatley.ClusterWheatleySoftDerivativesMix;
import etomica.virial.wheatley.ClusterWheatleySoftDerivativesMixBD;

import java.awt.*;
import java.util.*;
import java.util.List;

/**
 * Compute pure, binary, ternary and quaternary mixture virial coefficients using overlap sampling simulations
 * for some molecules using the TraPPE force fields.
 *
 */
public class VirialTraPPE {


    public static void main(String[] args) {
        VirialParam params = new VirialParam();
        VirialAlkane alkane = new VirialAlkane();
//        alkane.main(args);
        boolean isCommandline = args.length > 0;
        if (isCommandline) {
            ParseArgs.doParseArgs(params, args);
        } else {
            // Customize Interactive Parameters Here
            params.chemForm = new ChemForm[]{ChemForm.CH3OH};
            params.nPoints = 5; //B order
            params.nTypes = new int[]{5};
            params.nDer = 2;
            params.temperature = 1000;
            params.diagram = "808c";
            params.numSteps = 10000000;
            params.refFrac = -1;
            params.seed = null;
            params.dorefpref = false;
            params.doChainRef = true;

            params.BDtol = 1e-12;
        }

        // Import Params
        final ChemForm[] chemForm = params.chemForm;
        final int nPoints = params.nPoints;
        final int[] nTypes = params.nTypes;
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
        final double HSBn = doChainRef ? SpecialFunctions.factorial(nPoints) / 2 * Math.pow(vhs, nPoints - 1) : Standard.BHS(nPoints, sigmaHSRef);
        System.out.println("Chemform Length:" + chemForm.length);
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

        ClusterWheatleySoftDerivativesMixBD targetClusterBDRigid = null;
        SpeciesManager.Builder sb = SpeciesManager.builder();

        boolean anyPolar = false;
        MayerFunction[][] fAll = new MayerFunction[nTypes.length][nTypes.length];
        TraPPEParams[] TPList = new TraPPEParams[chemForm.length];
        for (int i = 0; i < TPList.length; i++) {
            TPList[i] = new TraPPEParams(space, chemForm[i]);
        }

        for (int i = 0; i < chemForm.length; i++) {
            sb.addSpecies(TPList[i].species);
        }
        SpeciesManager sm = sb.build();

        for (int i = 0; i < chemForm.length; i++) {
            TraPPEParams TPi = TPList[i];
            TPi.buildPotentials(sm);
            PotentialMoleculePair PGii = TPi.potentialGroup;

            P2PotentialGroupBuilder.ModelParams MPi = new P2PotentialGroupBuilder.ModelParams(TPi.atomTypes, TPi.sigma, TPi.epsilon, TPi.charge);
            fAll[i][i] = new MayerGeneral(PGii);

            anyPolar = (anyPolar || TPi.polar);

            for (int j = i + 1; j < chemForm.length; j++) {

                TraPPEParams TPj = TPList[j];

                P2PotentialGroupBuilder.ModelParams MPj = new P2PotentialGroupBuilder.ModelParams(TPj.atomTypes, TPj.sigma, TPj.epsilon, TPj.charge);

                PotentialMoleculePair PGij = P2PotentialGroupBuilder.P2PotentialGroupBuilder(space, sm, MPi, MPj);

                fAll[i][j] = fAll[j][i] = new MayerGeneral(PGij);

            }
        }

        //flex moves
        PotentialMasterBonding.FullBondingInfo bondingInfo = new PotentialMasterBonding.FullBondingInfo(sm);
        int nSpheres = TPList[0].species.getAtomTypes().size();
        boolean isFlex = TPList[0].isFlex && (params.diagram == null|| !params.diagram.equals("BC"));

        System.out.println("isFlex = " + isFlex);
        System.out.println("Diagram " + params.diagram);
        VirialDiagrams Diagrams = new VirialDiagrams(nPoints, false, isFlex);

        Diagrams.setDoReeHoover(false);
        ClusterAbstract targetCluster = Diagrams.makeVirialCluster(fAll[0][0]);
        ClusterSumShell[] targetDiagrams = new ClusterSumShell[0];
        int[] targetDiagramNumbers = new int[0];
        boolean[] diagramFlexCorrection = null;

        if (nSpheres > 2) {
            targetDiagrams = Diagrams.makeSingleVirialClusters((ClusterSum)targetCluster, null, fAll[0][0]);
            Set<Graph> singleGraphs = Diagrams.getMSMCGraphs(true, false);
            Map<Graph,Graph> cancelMap = Diagrams.getCancelMap();
            if(params.diagram != null && !params.diagram.equals("BC")) {
                int iGraph = 0;
                for (Graph g : singleGraphs) {
                    if(params.diagram.equals(g.getStore().toNumberString() + "c")) {
                        targetCluster = Diagrams.makeVirialCluster(g, fAll[0][0]);
                        targetDiagrams = new ClusterSumShell[0];
                        System.out.print(iGraph+" ("+g.coefficient()+") "+g.getStore().toNumberString()); // toNumberString: its corresponding number
                        Graph cancelGraph = cancelMap.get(g);
                        Set<Graph> gSplit = Diagrams.getSplitDisconnectedVirialGraphs(cancelGraph);

                        System.out.print(" - "+alkane.getSplitGraphString(gSplit, Diagrams, false));

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
                    targetDiagramNumbers[iGraph] = Integer.parseInt(g.getStore().toNumberString());

                    Graph cancelGraph = cancelMap.get(g);
                    if (cancelGraph != null) {
                        diagramFlexCorrection[iGraph] = true;
                        Set<Graph> gSplit = Diagrams.getSplitDisconnectedVirialGraphs(cancelGraph);

                        System.out.print(" - " + alkane.getSplitGraphString(gSplit, Diagrams, false));

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
                    System.out.println(g.coefficient()+" "+alkane.getSplitGraphString(gSplit, Diagrams, true));
                }
                System.out.println();
            }
        }

        targetCluster.setTemperature(temperature);
        for (int i=0; i<targetDiagrams.length; i++) {
            targetDiagrams[i].setTemperature(temperature);
        }


        //P3 bond angle
        if (TPList[0].theta_eq != null) {
            P3BondAngle[] p3 = new P3BondAngle[TPList[0].theta_eq.length]; //declaration, instatation

            List<int[]> triplets = new ArrayList<>();
            for (int i = 0; i < nSpheres - 2; i++) {
                triplets.add(new int[]{i, i + 1, i + 2});
            }
            for (int j = 0; j < TPList[0].theta_eq.length; j++) {
                p3[j] = new P3BondAngle(TPList[0].theta_eq[j], TPList[0].k_theta[j]);
                bondingInfo.setBondingPotentialTriplet(TPList[0].species, p3[j], Collections.singletonList(triplets.get(j)));

            }

        }
        //dihedral stuff
        P4BondTorsion[] p4 = null;
        if (TPList[0].a != null) {
            p4 = new P4BondTorsion[TPList[0].a.length];

            List<int[]> quads = new ArrayList<>();
            for (int i=0; i<nSpheres-3; i++) {
                quads.add(new int[]{i,i+1,i+2,i+3});
            }
            for (int i=0; i < TPList[0].a.length; i++) {
                p4[i] = new P4BondTorsion(space, TPList[0].a[i][0], TPList[0].a[i][1], TPList[0].a[i][2], TPList[0].a[i][3]);
                bondingInfo.setBondingPotentialQuad(TPList[0].species, p4[i], Collections.singletonList(quads.get(i)));

            }
            System.out.println(Arrays.toString(bondingInfo.bondedQuadPartners));
//            System.exit(0);

        }


        // Setting up Target Cluster for Rigid
        targetClusterRigid = new ClusterWheatleySoftDerivativesMix(nPoints, nTypes, fAll, BDtol, nDer);
        targetClusterRigid.setTemperature(temperature);

        // Setting BlockSize
        long blockSize = steps/numBlocks;
        int EqSubSteps = 1000;

        System.out.println(steps + " steps (" + numBlocks + " blocks of " + blockSize + ")");
        System.out.println("BD_Tol: " + BDtol + " BDAccFrac: " + BDAccFrac);

        // Setting up Flipping rigid, polar
        if(anyPolar && !isFlex && nPoints==2) {
            System.out.println("Performing Flipping");
            ((ClusterWheatleySoftDerivativesMix) targetClusterRigid).setTolerance(0);
            final int precision = -3*(int)Math.log10(BDtol);
            targetClusterBDRigid = new ClusterWheatleySoftDerivativesMixBD(nPoints,nTypes,fAll,precision,nDer);
            targetClusterBDRigid.setTemperature(temperature);
            ((ClusterWheatleySoftDerivativesMix) targetClusterRigid).setDoCaching(false);
            targetClusterBDRigid.setDoCaching(false);
            targetClusterRigid = new ClusterCoupledFlippedMultivalue(targetClusterRigid, targetClusterBDRigid, space, 20, nDer, BDtol);
        }
        //flipping for flexible polar B2
        else if(isFlex && nPoints==2 && anyPolar){
            targetCluster = new ClusterCoupledFlipped(targetCluster, space);

        }
        else if(anyPolar && isFlex && nPoints > 2 && params.diagram != null && !params.diagram.equals("BC")){
            int[][] flipPoints = Diagrams.getFlipPointsforDiagram(params.diagram);
            targetCluster = new ClusterCoupledFlippedPoints(targetCluster, space, flipPoints);

        }



        // Setting up Simulation
        SimulationVirialOverlap2 sim = null;

        if(!isFlex) {
            sim = new SimulationVirialOverlap2(space, sm, nTypes, temperature, refCluster, targetClusterRigid);
            if(seed!=null)sim.setRandom(new RandomMersenneTwister(seed));
            System.out.println("random seeds: "+ Arrays.toString(seed==null?sim.getRandomSeeds():seed));

            if (targetClusterRigid instanceof ClusterCoupledFlippedMultivalue) {
                ((ClusterCoupledFlippedMultivalue) targetClusterRigid).setBDAccFrac(BDAccFrac, sim.getRandom());
            } else {
                ((ClusterWheatleySoftDerivativesMix) targetClusterRigid).setBDAccFrac(BDAccFrac, sim.getRandom());
                ((ClusterWheatleySoftDerivativesMix) targetClusterRigid).setNumBDCheckBins(8);
            }
            // Adding derivative clusters to simulation
            ClusterMultiToSingle[] primes = new ClusterMultiToSingle[nDer];
            for(int m=0;m<primes.length;m++){
                primes[m] = new ClusterMultiToSingle(targetClusterRigid, m + 1);
            }
            sim.setExtraTargetClusters(primes);

        }
        else{
            sim = new SimulationVirialOverlap2(space, sm, new int[]{(nPoints+1)}, temperature, refCluster, targetCluster);
            sim.setExtraTargetClusters(targetDiagrams);
//            sim.setDoWiggle(nSpheres > 2);
            sim.setBondingInfo(bondingInfo);
            sim.setIntraPairPotentials(TPList[0].potentialGroup.getAtomPotentials());
            sim.setRandom(new RandomMersenneTwister(2));
        }
        // Initialize Simulation
        sim.init();

        // Set Position Definitions
        sim.box[0].setPositionDefinition(new MoleculePositionCOM(space));
        sim.box[1].setPositionDefinition(new MoleculePositionCOM(space));
//        sim.integrators[1].dodebug = true;
//        System.out.println(targetCluster.value(sim.box[1]));
//        System.exit(0);
        // Setting Chain Ref Moves
        IntArrayList[] bonding = new IntArrayList[3];
        bonding[0] = new IntArrayList(new int[]{1});
        for (int i=1; i<2; i++) {
            bonding[i] = new IntArrayList(new int[]{i-1,i+1});
        }
        bonding[2] = new IntArrayList(new int[]{1});

        int[] constraintMap = new int[nPoints+1];
        MCMoveClusterAngleBend mcMoveAngle = null;
        MCMoveClusterAngleBend mcMoveAngle1 = null;
        if (isFlex) {
            for (int i=0; i<nPoints; i++) {
                constraintMap[i] = i;
            }
            constraintMap[nPoints] = 0;
            ((MCMoveClusterMoleculeMulti)sim.mcMoveTranslate[1]).setConstraintMap(constraintMap);
            ((MCMoveClusterRotateMoleculeMulti)sim.mcMoveRotate[0]).setConstraintMap(constraintMap);

            ((MCMoveClusterRotateMoleculeMulti)sim.mcMoveRotate[1]).setConstraintMap(constraintMap);


        }
        if(TPList[0].theta_eq != null) {
            mcMoveAngle = new MCMoveClusterAngleBend(sim.getRandom(), sim.integrators[0].getPotentialCompute(), space);
            sim.integrators[0].getMoveManager().addMCMove(mcMoveAngle);
            mcMoveAngle1 = new MCMoveClusterAngleBend(sim.getRandom(), sim.integrators[1].getPotentialCompute(), space);
            sim.integrators[1].getMoveManager().addMCMove(mcMoveAngle1);

        }
        if (doChainRef) {
                sim.integrators[0].getMoveManager().removeMCMove(sim.mcMoveTranslate[0]);
//                sim.integrators[1].getMoveManager().removeMCMove(sim.mcMoveTranslate[1]);
//                sim.integrators[0].getMoveManager().removeMCMove(sim.mcMoveRotate[0]);

            MCMoveClusterMoleculeHSChain mcMoveHSC = new MCMoveClusterMoleculeHSChain(sim.getRandom(), sim.box[0], sigmaHSRef);
                if(isFlex) {
                    mcMoveHSC.setConstraintMap(constraintMap);

                }
                sim.integrators[0].getMoveManager().addMCMove(mcMoveHSC);
                sim.accumulators[0].setBlockSize(1);

        }


        // create the intramolecular potential here, add to it and add it to
        // the potential master if needed
        MCMoveClusterTorsionMulti[] torsionMoves = null;
        if (TPList[0].a != null) {

            torsionMoves = new MCMoveClusterTorsionMulti[TPList[0].a.length];
            for (int i=0; i<TPList[0].a.length; i++) {
                torsionMoves[i] = new MCMoveClusterTorsionMulti(sim.integrators[0].getPotentialCompute(), space, sim.getRandom(), p4[i], 40);
                torsionMoves[i].setBox(sim.box[0]);
                torsionMoves[i].setTemperature(temperature);
                sim.integrators[0].getMoveManager().addMCMove(torsionMoves[i]);
            }
//            torsionMoves[1] = new MCMoveClusterTorsionMulti(sim.integrators[1].getPotentialCompute(), space, sim.getRandom(), p4[0], 40);
//            torsionMoves[1].setBox(sim.box[1]);
//            torsionMoves[1].setTemperature(temperature);
//            sim.integrators[1].getMoveManager().addMCMove(torsionMoves[1]);
        }
        if(params.diagram != null && !params.diagram.equals("BC") && TPList[0].theta_eq != null) {

            ConfigurationClusterMoveMolecule move = new ConfigurationClusterMoveMolecule(space, sim.getRandom(), 5, new MCMoveBox[]{mcMoveAngle1});
            move.initializeCoordinates(sim.box[1]);
        }

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
            refFileName = "refpref_"+"_"+nPoints+"_"+tempString+"K";
        }

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
//        System.out.println(sim.integrators[0].getMoveManager().getMCMoves().get(2).getChi(150));


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
                System.out.println("Torsion move acceptance "+torsionMoves[i].getTracker().acceptanceRatio());}
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
        N2, O2, CO2, NH3, CH4, CH3OH, ethanol, propan1ol, propan2ol, isobutanol, C6H6, ethane, propane, ethaneUA, propaneUA, butaneUA, methaneUA, ethene, propene, butadiene, toluene
    }

    /**
     * Inner class for parameters
     */
    public static class VirialParam extends ParameterBase {
        // don't change these
        public ChemForm[] chemForm = {ChemForm.N2};
        public int nPoints = 2;
        public int[] nTypes = {2};
        public int nDer = 3;
        public double temperature = 400;
        public long numSteps = 1000000;

        public double refFrac = -1;
        public double sigmaHSRef = 5;
        public int[] seed = null;

        public boolean dorefpref = true;
        public boolean doChainRef = true;
        public String diagram = null;
        public double BDtol = 1e-12;

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
            else if (chemForm == ChemForm.CH4) {
                //TraPPE-EH
                //Atom in Compound
                //Avogadro 3D checked
                AtomType typeC = new AtomType(Carbon.INSTANCE);
                AtomType typeM = new AtomType(elementM);
                isFlex = false;
                atomTypes = new AtomType[]{typeC,typeM};

                //TraPPE Parameters
//                double bondLengthCM = 0.55; // Angstrom
                double sigmaC = 3.31; // Angstrom
                double epsilonC = Kelvin.UNIT.toSim(0.01);
                double qC = Electron.UNIT.toSim(0.0);
                double sigmaM = 3.31; // Angstrom
                double epsilonM = Kelvin.UNIT.toSim(15.30);
                double qM = Electron.UNIT.toSim(0.000);

                //Construct Arrays
                sigma = new double[] {sigmaC,sigmaM};
                epsilon = new double[] {epsilonC,epsilonM};
                charge = new double[]{qC, qM};

                //Get Coordinates
                Vector3D posC = new Vector3D(new double[]{0, 0, 0});
                Vector3D posM1 = new Vector3D(new double[]{0.0000000000, -0.4490790744,     -0.3175342264});
                Vector3D posM2 = new Vector3D(new double[]{0.0000000000,      0.4490790744,     -0.3175342264});
                Vector3D posM3 = new Vector3D(new double[]{-0.4490790744,      0.0000000000,      0.3175342264});
                Vector3D posM4 = new Vector3D(new double[]{0.4490790744,      0.0000000000,      0.3175342264  });

                //Set Geometry
                species = new SpeciesBuilder(space)
                        .addAtom(typeC, posC, "C")
                        .addAtom(typeM, posM1, "M1")
                        .addAtom(typeM, posM2, "M2")
                        .addAtom(typeM, posM3, "M3")
                        .addAtom(typeM, posM4, "M4")
                        .build();
            }
            else if (chemForm == ChemForm.CH3OH) {
                //TraPPE-UA
                //Atom in Compound
                //Avogadro 3D coordinates
                AtomType typeCH3 = new AtomType(Carbon.INSTANCE);
                AtomType typeO = new AtomType(Oxygen.INSTANCE);
                AtomType typeH = new AtomType(Hydrogen.INSTANCE);

                atomTypes = new AtomType[]{typeCH3,typeO, typeH};
                polar = true;
                isFlex = true;
                //TraPPE Parameters
//                double bondLengthCH3OH = 1.43; // Angstrom
//                double bondLengthOH = 0.945; //Angstrom
//                double thetaCH3OH = Degree.UNIT.toSim(18.5) ;
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


                //Construct Arrays
                sigma = new double[] {sigmaCH3,sigmaO, sigmaH};
                epsilon = new double[] {epsilonCH3,epsilonO, epsilonH};
                charge = new double[]{qCH3, qO, qH};
                theta_eq = new double[]{theta_CH3OH};
                k_theta = new double[]{k_thetaCH3OH};

                //Get Coordinates
                Vector3D posO = new Vector3D(new double[]{-0.6997355929,      0.0480707122,      0.0000000000 });
                Vector3D posCH3 = new Vector3D(new double[]{ 0.7278849876,     -0.0343882362,      0.0000000000});
                Vector3D posH = new Vector3D(new double[]{ -1.0507657087,     -0.8293133872,      0.0000000000   });
                System.out.println("pos0: "+ posO + ", posCH3: "+ posCH3 + ", posH: "+ posH);


                //Set Geometry
                species = new SpeciesBuilder(space)
                        .addAtom(typeCH3, posCH3, "CH3")
                        .addAtom(typeO, posO, "O")
                        .addAtom(typeH, posH, "H")
                        .build();
            }
            else if (chemForm == ChemForm.ethanol) {
                //TraPPE-UA
                //Atom in Compound
                //Avogadro 3d coordinates
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
                double thetaCH3OH = Degree.UNIT.toSim(18.5) ;
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
                sigma = new double[] {sigmaCH3,sigmaO, sigmaCH2, sigmaH};
                epsilon = new double[] {epsilonCH3,epsilonCH2, epsilonO, epsilonH};
                charge = new double[]{qCH3,qCH2, qO, qH};
                theta_eq = new double[]{theta_CCOH, theta_CH3OH};
                k_theta = new double[]{k_thetaCCOH, k_thetaCH3OH};
                a = new double[][]{{c0, c1, c2, c3}};
                double x3 = bondLengthCC + bondLengthCH3OH * Math.cos(theta_CCOH);
                double y3 = bondLengthCH3OH * Math.sin(theta_CCOH);
                double xH = x3 + bondLengthOH * Math.cos(theta_CH3OH);
                double yH = y3 + bondLengthOH * Math.sin(theta_CH3OH);
                //Get Coordinates
                Vector3D posCH2 = new Vector3D(new double[]{0.0343104921,     -0.5849517688,      0.0000000000 });
                Vector3D posCH3 = new Vector3D(new double[]{-1.2667321160 ,     0.2389949115,      0.0000000000});
                Vector3D posO = new Vector3D(new double[]{1.1583322479,      0.2990590299,      0.0000000000});

                Vector3D posH = new Vector3D(new double[]{1.9480255878,     -0.2199876240,      0.0000000000});

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
                //Avogadro 3D for dihedral combo (180, -60)
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
                sigma = new double[] {sigmaCH3,sigmaO, sigmaCH2_2, sigmaCH2_3, sigmaH};
                epsilon = new double[] {epsilonCH3,epsilonCH2_2, epsilonCH2_3, epsilonO, epsilonH};
                charge = new double[]{qCH3,qCH2_2, qCH2_3, qO, qH};
                theta_eq = new double[]{theta_CCC, theta_CCOH, theta_CH3OH};
                k_theta = new double[]{k_thetaCCC, k_thetaCCOH, k_thetaCH3OH};
                a = new double[][]{{c00, c01, c02, c03}, {c10, c11, c12, c13}};



                double x3 = bondLengthCC + bondLengthCC * Math.cos(theta_CCC);
                double y3 = bondLengthCC * Math.sin(theta_CCC);
                double xO = x3 + bondLengthCH3OH * Math.cos(theta_CCOH);
                double yO = y3 + bondLengthCH3OH * Math.sin(theta_CCOH);
                double xH = xO + bondLengthOH * Math.cos(theta_CH3OH);
                double yH = yO + bondLengthOH * Math.sin(theta_CH3OH);
                //Get Coordinates
                Vector3D posCH3 = new Vector3D(new double[]{-0.5506818169,      1.4780955839,     -1.2182943870});
                Vector3D posCH2_2 = new Vector3D(new double[]{-0.5697691905,      0.5363897718,      0.0000757684});
                Vector3D posCH2_3 = new Vector3D(new double[]{0.5592739627,     -0.5109218115,      0.0000000000});

                Vector3D posO = new Vector3D(new double[]{0.4503423313,     -1.3293023461,      1.1675989038});

                Vector3D posH = new Vector3D(new double[]{0.5241155912,     -0.7755909625,      1.9298234179});

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
                //Avogadro 3D structure
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
                sigma = new double[] {sigmaCH3,sigmaO, sigmaCH, sigmaH};
                epsilon = new double[] {epsilonCH3,epsilonCH, epsilonO, epsilonH};
                charge = new double[]{qCH3,qCH, qO, qH};
                theta_eq = new double[]{theta_CCC};
                k_theta = new double[]{k_thetaCCC};
                a = new double[][]{{c00, c01, c02, c03}}; //1-2-4-5, not 1-2-3-4, fix



                double x3 = bondLengthCC + bondLengthCC * Math.cos(theta_CCC);
                double y3 = bondLengthCC * Math.sin(theta_CCC);

                double xO = bondLengthCC + bondLengthCH3OH * Math.sin(theta_CCOH);
                double yO = bondLengthCH3OH * Math.cos(theta_CCOH);

                double xH = xO + bondLengthOH * Math.cos(theta_CH3OH);
                double yH = yO + bondLengthOH * Math.sin(theta_CH3OH);


                //Get Coordinates
                Vector3D posCH3 = new Vector3D(new double[]{-0.0134659706,     -0.0495590344,     -1.9255282468});
                Vector3D posCH = new Vector3D(new double[]{0.0052304145,     -0.0136988872,     -0.3860593465});
                Vector3D posCH3_3 = new Vector3D(new double[]{-0.5947251828,      1.2919716877,      0.1679060235});

                Vector3D posO = new Vector3D(new double[]{-0.7340205938,     -1.1293415899,      0.1176760074});

                Vector3D posH = new Vector3D(new double[]{-1.6329601728,     -1.0438205468,     -0.1609254724});

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
                //built own compound via ChemDoodle
                //Avogadro 3D compound dihedrals (180, 57)
                //https://pubs.acs.org/doi/epdf/10.1021/ja5011288?ref=article_openPDF
                //Atom in Compound
                AtomType typeCH3 = new AtomType(Carbon.INSTANCE);
                AtomType typeCH = new AtomType(Carbon.INSTANCE);
                AtomType typeCH2 = new AtomType(Carbon.INSTANCE);

                AtomType typeO = new AtomType(Oxygen.INSTANCE);
                AtomType typeH = new AtomType(Hydrogen.INSTANCE);

                atomTypes = new AtomType[]{typeCH3, typeCH, typeCH2, typeO, typeH};
                polar = true;
                isFlex = true;
                //TraPPE Parameters
//                double bondLengthCC= 1.54; //Angstrom
//                double bondLengthCH3OH = 1.43; // Angstrom
//                double bondLengthOH = 0.945; //Angstrom
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


                //Construct Arrays
                sigma = new double[] {sigmaCH3,sigmaO, sigmaCH, sigmaCH2, sigmaH};
                epsilon = new double[] {epsilonCH3,epsilonCH, epsilonO, epsilonCH2, epsilonH};
                charge = new double[]{qCH3,qCH, qO, qH, qCH2};
                theta_eq = new double[]{theta_CCC};
                k_theta = new double[]{k_thetaCCC};
                a = new double[][]{{c00, c01, c02, c03}}; //1-2-4-5, not 1-2-3-4, fix



//                double x3 = bondLengthCC + bondLengthCC * Math.cos(theta_CCC);
//                double y3 = bondLengthCC * Math.sin(theta_CCC);
//                double xO = bondLengthCC + bondLengthCH3OH * Math.cos(theta_CCOH);
//                double yO = bondLengthCH3OH * Math.sin(theta_CCOH);
//                double xH = xO + bondLengthOH * Math.cos(theta_CH3OH);
//                double yH = yO + bondLengthOH * Math.sin(theta_CH3OH);
                //Get Coordinates
                Vector3D posC1 = new Vector3D(new double[]{-0.7325836362,      1.4388713493,      0.0029000279});
                Vector3D posCH = new Vector3D(new double[]{-0.4954281363,     -0.0316338865,      0.3940158934});
                Vector3D posCH2 = new Vector3D(new double[]{0.7486712378,     -0.6191587529,     -0.2978159349});
                Vector3D posC4 = new Vector3D(new double[]{-1.7147587460,     -0.9169711829,      0.0761940314});

                Vector3D posO = new Vector3D(new double[]{0.5135119555,     -0.7088722652,     -1.7054918911});
                Vector3D posH = new Vector3D(new double[]{1.2946829447,     -1.0428618814,     -2.1193128934});

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



            else if (chemForm == ChemForm.C6H6) {
                //TraPPE-EH
                //Atom in Compound
                //planar
                AtomType typeC = new AtomType(Carbon.INSTANCE);
                AtomType typeH = new AtomType(Hydrogen.INSTANCE);

                atomTypes = new AtomType[]{typeC, typeH};
                isFlex = false;
                //TraPPE Parameters
                double bondLengthCC = 1.392; // Angstrom
                double bondLengthCH = 1.08; //Angstrom
                double sumbond = bondLengthCC + bondLengthCH;
                double theta = Degree.UNIT.toSim(60);
                double sigmaC = 3.60; // Angstrom
                double epsilonC = Kelvin.UNIT.toSim(30.70);
                double qC = Electron.UNIT.toSim(-0.095);
                double sigmaH = 2.36; // Angstrom
                double epsilonH = Kelvin.UNIT.toSim(25.45);
                double qH = Electron.UNIT.toSim(0.095);


                //Construct Arrays
                sigma = new double[] {sigmaC,sigmaH};
                epsilon = new double[] {epsilonC,epsilonH};
                charge = new double[]{qC, qH};

                //Get Coordinates
                Vector3D posC1 = new Vector3D(new double[]{bondLengthCC,0,  0});
                Vector3D posC2 = new Vector3D(new double[]{bondLengthCC * Math.cos(theta), bondLengthCC * Math.sin(theta), 0});
                Vector3D posC3 = new Vector3D(new double[]{-bondLengthCC * Math.cos(theta), bondLengthCC * Math.sin(theta), 0});
                Vector3D posC4 = new Vector3D(new double[]{-bondLengthCC,0, 0});
                Vector3D posC5 = new Vector3D(new double[]{-bondLengthCC * Math.cos(theta), -bondLengthCC * Math.sin(theta), 0});
                Vector3D posC6 = new Vector3D(new double[]{bondLengthCC * Math.cos(theta), -bondLengthCC * Math.sin(theta), 0});

                Vector3D posH1 = new Vector3D(new double[]{sumbond,0,  0});
                Vector3D posH2 = new Vector3D(new double[]{sumbond * Math.cos(theta), sumbond * Math.sin(theta), 0});
                Vector3D posH3 = new Vector3D(new double[]{-sumbond * Math.cos(theta), sumbond * Math.sin(theta), 0});
                Vector3D posH4 = new Vector3D(new double[]{-sumbond,0, 0});
                Vector3D posH5 = new Vector3D(new double[]{-sumbond * Math.cos(theta), -sumbond * Math.sin(theta), 0});
                Vector3D posH6 = new Vector3D(new double[]{sumbond * Math.cos(theta), -sumbond * Math.sin(theta), 0});
                System.out.println("Carbon Positions: " + Arrays.toString(new Vector3D[]{posC1, posC2, posC3, posC4, posC5, posC6}));
                System.out.println("Hydrogen Positions: " + Arrays.toString(new Vector3D[]{posH1, posH2, posH3, posH4, posH5, posH6}));



                //Set Geometry
                species = new SpeciesBuilder(space)
                        .addAtom(typeC, posC1, "C1")
                        .addAtom(typeC, posC2, "C2")
                        .addAtom(typeC, posC3, "C3")
                        .addAtom(typeC, posC4, "C4")
                        .addAtom(typeC, posC5, "C5")
                        .addAtom(typeC, posC6, "C6")

                        .addAtom(typeH, posH1, "H1")
                        .addAtom(typeH, posH2, "H2")
                        .addAtom(typeH, posH3, "H3")
                        .addAtom(typeH, posH4, "H4")
                        .addAtom(typeH, posH5, "H5")
                        .addAtom(typeH, posH6, "H6")

                        .build();
            }
            else if (chemForm == ChemForm.toluene) {
                    //TraPPE-UA
                //DO THIS
                //Atom in Compound
                //planar
                AtomType typeC = new AtomType(Carbon.INSTANCE);
                AtomType typeH = new AtomType(Hydrogen.INSTANCE);

                atomTypes = new AtomType[]{typeC, typeH};
                isFlex = false;
                //TraPPE Parameters
                double bondLengthCC = 1.392; // Angstrom
                double bondLengthCH = 1.08; //Angstrom
                double sumbond = bondLengthCC + bondLengthCH;
                double theta = Degree.UNIT.toSim(60);
                double sigmaC = 3.60; // Angstrom
                double epsilonC = Kelvin.UNIT.toSim(30.70);
                double qC = Electron.UNIT.toSim(-0.095);
                double sigmaH = 2.36; // Angstrom
                double epsilonH = Kelvin.UNIT.toSim(25.45);
                double qH = Electron.UNIT.toSim(0.095);


                //Construct Arrays
                sigma = new double[] {sigmaC,sigmaH};
                epsilon = new double[] {epsilonC,epsilonH};
                charge = new double[]{qC, qH};

                //Get Coordinates
                Vector3D posC1 = new Vector3D(new double[]{bondLengthCC,0,  0});
                Vector3D posC2 = new Vector3D(new double[]{bondLengthCC * Math.cos(theta), bondLengthCC * Math.sin(theta), 0});
                Vector3D posC3 = new Vector3D(new double[]{-bondLengthCC * Math.cos(theta), bondLengthCC * Math.sin(theta), 0});
                Vector3D posC4 = new Vector3D(new double[]{-bondLengthCC,0, 0});
                Vector3D posC5 = new Vector3D(new double[]{-bondLengthCC * Math.cos(theta), -bondLengthCC * Math.sin(theta), 0});
                Vector3D posC6 = new Vector3D(new double[]{bondLengthCC * Math.cos(theta), -bondLengthCC * Math.sin(theta), 0});

                Vector3D posH1 = new Vector3D(new double[]{sumbond,0,  0});
                Vector3D posH2 = new Vector3D(new double[]{sumbond * Math.cos(theta), sumbond * Math.sin(theta), 0});
                Vector3D posH3 = new Vector3D(new double[]{-sumbond * Math.cos(theta), sumbond * Math.sin(theta), 0});
                Vector3D posH4 = new Vector3D(new double[]{-sumbond,0, 0});
                Vector3D posH5 = new Vector3D(new double[]{-sumbond * Math.cos(theta), -sumbond * Math.sin(theta), 0});
                Vector3D posH6 = new Vector3D(new double[]{sumbond * Math.cos(theta), -sumbond * Math.sin(theta), 0});
                System.out.println("Carbon Positions: " + Arrays.toString(new Vector3D[]{posC1, posC2, posC3, posC4, posC5, posC6}));
                System.out.println("Hydrogen Positions: " + Arrays.toString(new Vector3D[]{posH1, posH2, posH3, posH4, posH5, posH6}));



                //Set Geometry
                species = new SpeciesBuilder(space)
                        .addAtom(typeC, posC1, "C1")
                        .addAtom(typeC, posC2, "C2")
                        .addAtom(typeC, posC3, "C3")
                        .addAtom(typeC, posC4, "C4")
                        .addAtom(typeC, posC5, "C5")
                        .addAtom(typeC, posC6, "C6")

                        .addAtom(typeH, posH1, "H1")
                        .addAtom(typeH, posH2, "H2")
                        .addAtom(typeH, posH3, "H3")
                        .addAtom(typeH, posH4, "H4")
                        .addAtom(typeH, posH5, "H5")
                        .addAtom(typeH, posH6, "H6")

                        .build();
            }

            else if (chemForm == ChemForm.ethane) {
                //TraPPE-EH
                //Atom in Compound
                //avogadro
                AtomType typeC = new AtomType(Carbon.INSTANCE);
                AtomType typeM = new AtomType(elementM);

                atomTypes = new AtomType[]{typeC,typeM};
                isFlex = false;
                //TraPPE Parameters
                double bondLengthCM = 0.55; // Angstrom
                double bondLengthCC = 1.5350; //Angstrom
                double thetaCCH = Degree.UNIT.toSim(110.70);
                double a120 = Degree.UNIT.toSim(120);
                double a240 = Degree.UNIT.toSim(240);
                double sigmaC = 3.30; // Angstrom
                double epsilonC = Kelvin.UNIT.toSim(4.00);
                double qC = Electron.UNIT.toSim(0.0);
                double sigmaM = 3.31; // Angstrom
                double epsilonM = Kelvin.UNIT.toSim(15.30);
                double qM = Electron.UNIT.toSim(0.000);
                double a01 = Kelvin.UNIT.toSim(717);
                //Construct Arrays
                sigma = new double[] {sigmaC,sigmaM};
                epsilon = new double[] {epsilonC,epsilonM};
                charge = new double[]{qC, qM};
                a = new double[][]{{0, a01, 0, 0}};


                //Get Coordinates
                Vector3D posC1 = new Vector3D(new double[]{0.7674410120,     -0.0579174161  ,   -0.0637971165});
                Vector3D posC2 = new Vector3D(new double[]{-0.7632720540,      0.0567215202,     -0.0644482175});
                Vector3D posM1 = new Vector3D(new double[]{0.9474423656,     -0.3387486927,      0.3735051894});
                Vector3D posM2 = new Vector3D(new double[]{0.9394900365 ,    -0.3210710943,     -0.5150722819});
                Vector3D posM3 = new Vector3D(new double[]{0.9997317223,      0.4406138058,     -0.0610318522});
                Vector3D posM4 = new Vector3D(new double[]{-0.9436953068,      0.3368882196,      0.3731064310});
                Vector3D posM5 = new Vector3D(new double[]{-0.9348859373 ,     0.3205606343,     -0.5154887887});
                Vector3D posM6 = new Vector3D(new double[]{-0.9955649581,     -0.4418131574,     -0.0626644104});

                //Set Geometry
                species = new SpeciesBuilder(space)
                        .addAtom(typeC, posC1, "C1")
                        .addAtom(typeC, posC2, "C2")
                        .addAtom(typeM, posM1, "M1")
                        .addAtom(typeM, posM2, "M2")
                        .addAtom(typeM, posM3, "M3")
                        .addAtom(typeM, posM4, "M4")
                        .addAtom(typeM, posM5, "M5")
                        .addAtom(typeM, posM6, "M6")

                        .build();
            }
            else if (chemForm == ChemForm.propane) {
                //TraPPE-EH
                //Atom in Compound
                //Avogadro
                AtomType typeC = new AtomType(Carbon.INSTANCE);
                AtomType typemidC = new AtomType(Carbon.INSTANCE);
                AtomType typeM = new AtomType(elementM);

                atomTypes = new AtomType[]{typeC, typemidC, typeM};
                isFlex = true;
                //TraPPE Parameters
                double bondLengthCM = 0.55; // Angstrom
                double bondLengthCC = 1.5350; //Angstrom
                double alpha = Degree.UNIT.toSim(53.90);
                double thetaCCH = Degree.UNIT.toSim(110.70);
                double thetaCCC = Degree.UNIT.toSim(112.70);
                double a120 = Degree.UNIT.toSim(120);
                double a240 = Degree.UNIT.toSim(240);
                double sigmaC = 3.30; // Angstrom
                double epsilonC = Kelvin.UNIT.toSim(4.00);
                double qC = Electron.UNIT.toSim(0.0);
                double sigmamidC = 3.65; // Angstrom
                double epsilonmidC = Kelvin.UNIT.toSim(5.00);
                double qmidC = Electron.UNIT.toSim(0.0);
                double sigmaM = 3.31; // Angstrom
                double epsilonM = Kelvin.UNIT.toSim(15.30);
                double qM = Electron.UNIT.toSim(0.000);
                double k = Kelvin.UNIT.toSim(58765);
                double a01 = Kelvin.UNIT.toSim(854);
                //Construct Arrays
                sigma = new double[] {sigmaC,sigmamidC, sigmaM};
                epsilon = new double[] {epsilonC,epsilonmidC, epsilonM};
                charge = new double[]{qC, qmidC, qM};
                k_theta = new double[]{k};
                theta_eq = new double[]{thetaCCC};
                a = new double[][]{{0, a01, 0, 0}}; //fix, 3-1-2-4


                double x3 = bondLengthCC - bondLengthCC * Math.cos(thetaCCC);
                double y3 = bondLengthCC * Math.sin(thetaCCC);
                //Get Coordinates
                Vector3D posC1 = new Vector3D(new double[]{1.2520730479,     -0.2644906931,     -0.1533486375});
                Vector3D posC2 = new Vector3D(new double[]{0.0099234689,      0.6324679426,     -0.0597824612});
                Vector3D posC3 = new Vector3D(new double[]{-1.2958658885 ,    -0.1704726300,      0.0203773238});
                Vector3D posM1 = new Vector3D(new double[]{1.2785732026,     -0.6058207170,      0.2771064201});
                Vector3D posM2 = new Vector3D(new double[]{1.2323681877,     -0.5707199860,     -0.6097862988});
                Vector3D posM3 = new Vector3D(new double[]{1.7030456657 ,     0.0491621611,     -0.1806540407});
                Vector3D posM4 = new Vector3D(new double[]{0.0613163156,      0.9576922390,      0.3807716584});
                Vector3D posM5 = new Vector3D(new double[]{0.0003552066,      0.9539846380,     -0.5059163682});
                Vector3D posM6 = new Vector3D(new double[]{-1.3069583182,     -0.4786267642,      0.4758089735});
                Vector3D posM7 = new Vector3D(new double[]{-1.7235437146 ,     0.1747612163,      0.0405086015});
                Vector3D posM8 = new Vector3D(new double[]{-1.3538001794,     -0.4839632176,     -0.4278048595});

                //Set Geometry
                species = new SpeciesBuilder(space)
                        .addAtom(typeC, posC1, "C1")
                        .addAtom(typemidC, posC2, "C2")
                        .addAtom(typeC, posC3, "C3")
                        .addAtom(typeM, posM1, "M1")
                        .addAtom(typeM, posM2, "M2")
                        .addAtom(typeM, posM3, "M3")
                        .addAtom(typeM, posM4, "M4")
                        .addAtom(typeM, posM5, "M5")
                        .addAtom(typeM, posM6, "M6")
                        .addAtom(typeM, posM7, "M7")
                        .addAtom(typeM, posM8, "M8")
                        .build();
                System.out.println("Carbon Positions: " + Arrays.toString(new Vector3D[]{posC1, posC2, posC3}));
                System.out.println("Hydrogen Positions: " + Arrays.toString(new Vector3D[]{posM1, posM2, posM3, posM4, posM5, posM6, posM7, posM8}));
            }
            else if (chemForm == ChemForm.ethaneUA) {
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
            else if (chemForm == ChemForm.propaneUA) {
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
                //Construct Arrays
                sigma = new double[] {sigmaCH2, sigmaCH3};
                epsilon = new double[] {epsilonCH2, epsilonCH3};
                charge = new double[]{qCH2, qCH3};
                k_theta = new double[]{kCCC};
                theta_eq = new double[]{thetaCCH};
                //Get Coordinates
                Vector3D posC1 = new Vector3D(new double[]{1.2655314529,     -0.2810000000,      0.0000000000});
                Vector3D posC2 = new Vector3D(new double[]{0.0081657264,      0.6081745779,      0.0000000000});
                Vector3D posC3 = new Vector3D(new double[]{-1.3155523850,     -0.1788262871,      0.0000000000});

                //Set Geometry
                species = new SpeciesBuilder(space)
                        .addAtom(typeCH3, posC1, "C1")
                        .addAtom(typeCH2, posC2, "C2")
                        .addAtom(typeCH3, posC3, "C3")
                        .build();
                System.out.println("Positions:"+ posC1+ ',' +posC2+ ','+ posC3);

            }
            else if (chemForm == ChemForm.butaneUA) {
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
                sigma = new double[] {sigmaCH2, sigmaCH3};
                epsilon = new double[] {epsilonCH2, epsilonCH3};
                charge = new double[]{qCH2, qCH3};
                k_theta = new double[]{k, k};
                theta_eq = new double[]{thetaCCH, thetaCCH};
                a = new double[][]{{a00, a01, a02, a03}};

                double x3 = bondLengthCHxCHy - bondLengthCHxCHy * Math.cos(thetaCCH);
                double y3 = bondLengthCHxCHy * Math.sin(thetaCCH);
                //Get Coordinates


                Vector3D posC1 = new Vector3D(new double[]{1.9453970388,     -0.1390603157,      0.0000000000});
                Vector3D posC2 = new Vector3D(new double[]{0.5623293516,      0.5382311818,      0.0000000000});
                Vector3D posC3 = new Vector3D(new double[]{-0.6189515277,     -0.4497847514,      0.0000000000});
                Vector3D posC4 = new Vector3D(new double[]{-2.0020192149,      0.2275067461,      0.0000000000});

                System.out.println("Carbon Positions: " + Arrays.toString(new Vector3D[]{posC1, posC2, posC3, posC4}));

                //Set Geometry
                species = new SpeciesBuilder(space)
                        .addAtom(typeCH3, posC1, "C1")
                        .addAtom(typeCH2, posC2, "C2")
                        .addAtom(typeCH2, posC3, "C3")
                        .addAtom(typeCH3, posC4, "C4")

                        .build();
            }
            else if (chemForm == ChemForm.methaneUA) {
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
                //Avogadro
                AtomType typeCH3 = new AtomType(Carbon.INSTANCE);
                AtomType typeCH = new AtomType(Carbon.INSTANCE);
                AtomType typeCH2 = new AtomType(Carbon.INSTANCE);

                atomTypes = new AtomType[]{typeCH3, typeCH, typeCH2};
                isFlex = true;
                //TraPPE Parameters
                double bondLengthCHxCHy = 1.54; // Angstrom
                double bondLengthCHxChy_double = 1.33; //Angstrom
                double thetaCCC = Degree.UNIT.toSim(29.70);//thetaeq
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
                double x3 = bondLengthCHxChy_double * Math.sin(thetaCCC) + bondLengthCHxCHy;
                //Get Coordinates

                Vector3D posCH3 = new Vector3D(new double[]{1.2546345088,     -0.2478094308,     -0.0004026143});
                Vector3D posCH = new Vector3D(new double[]{0.1513145765,      0.4948780916,      0.0002065621});
                Vector3D posCH2 = new Vector3D(new double[]{-1.2286312755,     -0.1887514487,     -0.0001051992});

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
                //Avogardo planar molecule
                AtomType typeCH = new AtomType(Carbon.INSTANCE);
                AtomType typeCH2 = new AtomType(Carbon.INSTANCE);

                atomTypes = new AtomType[]{typeCH, typeCH2};
                isFlex = true;
                //TraPPE Parameters
                double bondLengthCHxCHy = 1.54; // Angstrom
                double bondLengthCHxChy_double = 1.33; //Angstrom
                double thetaCCC = Degree.UNIT.toSim(29.70);

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
                double x3 = bondLengthCHxCHy * Math.sin(thetaCCC) + bondLengthCHxChy_double;
                double y3 = bondLengthCHxCHy * Math.cos(thetaCCC);

                //Get Coordinates
                Vector3D posC1 = new Vector3D(new double[]{1.8660853064,     -0.0921588646,     -0.0001862610});
                Vector3D posC2 = new Vector3D(new double[]{0.6579550041,      0.4640075516,      0.0001137386});
                Vector3D posC3 = new Vector3D(new double[]{-0.5945193236,     -0.4320437626,     -0.0002041153});
                Vector3D posC4 = new Vector3D(new double[]{-1.8026496259,      0.1241226536,      0.0000958844});

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
            else if (chemForm == ChemForm.C6H6) {
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
            else if (chemForm == ChemForm.ethaneUA) {
                P2PotentialGroupBuilder.ModelParams modelParams = new P2PotentialGroupBuilder.ModelParams(atomTypes, sigma, epsilon, charge);
                potentialGroup = P2PotentialGroupBuilder.P2PotentialGroupBuilder(space, sm, modelParams, null);
            }
            else if (chemForm == ChemForm.propaneUA) {
                P2PotentialGroupBuilder.ModelParams modelParams = new P2PotentialGroupBuilder.ModelParams(atomTypes, sigma, epsilon, charge);
                potentialGroup = P2PotentialGroupBuilder.P2PotentialGroupBuilder(space, sm, modelParams, null);
            }
            else if (chemForm == ChemForm.butaneUA) {
                P2PotentialGroupBuilder.ModelParams modelParams = new P2PotentialGroupBuilder.ModelParams(atomTypes, sigma, epsilon, charge);
                potentialGroup = P2PotentialGroupBuilder.P2PotentialGroupBuilder(space, sm, modelParams, null);
            }
            else if (chemForm == ChemForm.methaneUA) {
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



        }

    }

}
