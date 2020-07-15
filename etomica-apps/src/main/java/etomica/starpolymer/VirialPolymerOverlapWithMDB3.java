/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.starpolymer;

import etomica.action.IAction;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.chem.elements.ElementSimple;
import etomica.data.IData;
import etomica.data.types.DataGroup;
import etomica.graph.model.Graph;
import etomica.graph.operations.DeleteEdge;
import etomica.graph.operations.DeleteEdgeParameters;
import etomica.graphics.*;
import etomica.integrator.mcmove.MCMoveStepTracker;
import etomica.math.SpecialFunctions;
import etomica.molecule.IMoleculeList;
import etomica.molecule.MoleculePositionCOM;
import etomica.potential.IPotential;
import etomica.potential.P2HardSphere;
import etomica.potential.PotentialGroup;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.virial.*;
import etomica.virial.cluster.VirialDiagrams;
import etomica.virial.simulations.SimulationVirialOverlap2;

import java.awt.*;
import java.util.Map;
import java.util.Set;

/**
 * Calculation for virial coefficients of hard spheres
 * Edited for calculating protein 1HRB
 * <p>
 * CHANGES: change for star-polymers. May 7, 2018
 */
public class VirialPolymerOverlapWithMDB3 {

    public static void main(String[] args) {

        VirialHSParam params = new VirialHSParam();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        } else {
            params.nPoints = 3;
            params.numSteps = 5000000L;
            params.ref = VirialHSParam.CHAINS;
            params.fv = 3;
            params.lv = 16;
        }

        final int nPoints = params.nPoints;
        long steps = params.numSteps;
        final int ref = params.ref;
        final int f = params.fv;
        final int l = params.lv;
        final double temperature = 1.0;

        Space space = Space3D.getInstance();
        SpeciesPolymerMono species = new SpeciesPolymerMono(space, new ElementSimple("A"), f, l);

        final double sigmaTranslate = 30.0;
        final double refDiameter = 11.0;
        double vhs = (4.0 / 3.0) * Math.PI * refDiameter * refDiameter * refDiameter;

        System.out.println("HS singly-connected sampling B" + nPoints);
        System.out.println("Computing for Star-polymer f = " + f + " l = " + l);

        final P2HardSphere potentialHS = new P2HardSphere(space);
        final PotentialGroup potentialGroup = new PotentialGroup(2, space);
        potentialGroup.addPotential(potentialHS, new AtomType[]{species.getAtomType(0), species.getAtomType(0)});

        MayerGeneral fRef = new MayerGeneral(potentialGroup);
        MayerFunction fRefPos = new MayerFunction() {
            public void setBox(Box box) {
            }

            public IPotential getPotential() {
                return null;
            }

            public double f(IMoleculeList pair, double r2, double beta) {
                return r2 < refDiameter * refDiameter ? 1 : 0;
            }
        };

        ClusterAbstract refCluster = null;
        long numDiagrams;
        double ri = 0;
        if (ref == VirialHSParam.CHAINS) {
            System.out.println("using a chain reference");
            refCluster = new ClusterChainHS(nPoints, fRefPos);
            numDiagrams = ((ClusterChainHS) refCluster).numDiagrams();
            ri = numDiagrams * Math.pow(vhs, nPoints - 1);
        }

        // (nPoints-1)! is simply not included by ClusterWheatley, so do that here.
        final double refIntegral = ri / SpecialFunctions.factorial(nPoints);
        System.out.println("reference integral: " + refIntegral);
//        refCluster.setTemperature(temperature);
        System.out.println(steps + " steps");

//        ClusterAbstract targetCluster = new ClusterWheatleyHS(nPoints, fRef);
//        targetCluster.setTemperature(temperature);

        // flexibility contribution starts here
        boolean flex = nPoints > 2;
        ClusterSumShell[] targetDiagrams = new ClusterSumShell[0];
        int[] targetDiagramNumbers = new int[0];
        boolean[] diagramFlexCorrection = null;

        MayerGeneral fTarget = new MayerGeneral(potentialGroup);
        VirialDiagrams polymerDiagrams = new VirialDiagrams(nPoints, false, flex);
        polymerDiagrams.setDoReeHoover(true);
        ClusterSum targetClusterFlex = polymerDiagrams.makeVirialCluster(fTarget);

//            VirialDiagrams rigidDiagrams = new VirialDiagrams(nPoints, false, false);
//            ClusterSum refCluster = rigidDiagrams.makeVirialCluster(fRef);
//            double refIntegral = HSB[nPoints];

        targetDiagrams = polymerDiagrams.makeSingleVirialClusters(targetClusterFlex, null, fTarget);
        targetDiagramNumbers = new int[targetDiagrams.length];
        System.out.println("individual clusters:");
        Set<Graph> singleGraphs = polymerDiagrams.getMSMCGraphs(true, false);
        Map<Graph, Graph> cancelMap = polymerDiagrams.getCancelMap();
        int iGraph = 0;
        diagramFlexCorrection = new boolean[targetDiagrams.length];
        for (Graph g : singleGraphs) {
            System.out.print(iGraph + " (" + g.coefficient() + ") " + g.getStore().toNumberString()); // toNumberString: its corresponding number
            targetDiagramNumbers[iGraph] = Integer.parseInt(g.getStore().toNumberString());

            Graph cancelGraph = cancelMap.get(g);
            if (cancelGraph != null) {
                diagramFlexCorrection[iGraph] = true;
                Set<Graph> gSplit = polymerDiagrams.getSplitDisconnectedVirialGraphs(cancelGraph);

                System.out.print(" - " + getSplitGraphString(gSplit, polymerDiagrams, false));

            }
            System.out.println();
            iGraph++;
        }
        System.out.println();
        Set<Graph> disconnectedGraphs = polymerDiagrams.getExtraDisconnectedVirialGraphs();
        if (disconnectedGraphs.size() > 0) {
            System.out.println("extra clusters:");

            for (Graph g : disconnectedGraphs) {
                Set<Graph> gSplit = polymerDiagrams.getSplitDisconnectedVirialGraphs(g);
                System.out.println(g.coefficient() + " " + getSplitGraphString(gSplit, polymerDiagrams, true));
            }
            System.out.println();
        }

        targetClusterFlex.setTemperature(temperature);
        for (ClusterSumShell targetDiagram : targetDiagrams) {
            targetDiagram.setTemperature(temperature);
        }

        ClusterWeight[] sampleClusters = new ClusterWeight[]{ClusterWeightAbs.makeWeightCluster(refCluster), ClusterWeightAbs.makeWeightCluster(targetClusterFlex)};

        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space, new ISpecies[]{species},
                new int[]{flex ? nPoints + 1 : nPoints}, temperature, new ClusterAbstract[]{refCluster, targetClusterFlex}, targetDiagrams, sampleClusters, false);

        int[] constraintMap = new int[nPoints + 1];
        for (int i = 0; i < nPoints; i++) {
            constraintMap[i] = i;
        }
        constraintMap[nPoints] = 0;
        ((MCMoveClusterMoleculeMulti) sim.mcMoveTranslate[0]).setConstraintMap(constraintMap);
        ((MCMoveClusterMoleculeMulti) sim.mcMoveTranslate[1]).setConstraintMap(constraintMap);
        ((MCMoveClusterRotateMoleculeMulti) sim.mcMoveRotate[0]).setConstraintMap(constraintMap);
        ((MCMoveClusterRotateMoleculeMulti) sim.mcMoveRotate[1]).setConstraintMap(constraintMap);

        sim.integratorOS.setRefStepFraction(-1);
        sim.integratorOS.setAdjustStepFraction(true);

//        if (ref == VirialHSParam.CHAINS && flex) {
//            sim.integrators[0].getMoveManager().removeMCMove(sim.mcMoveTranslate[0]);
//            sim.integrators[0].getMoveManager().addMCMove(new MCMoveClusterMoleculeHSChainB3(sim.getRandom(), space, refDiameter));
//            sim.accumulators[0].setBlockSize(1);
//        }

        ((MCMoveClusterRotateMoleculeMulti) sim.mcMoveRotate[0]).setPositionDefinition(new MoleculePositionCOM(space));
//        sim.integrators[1].getMoveManager().removeMCMove(sim.mcMoveRotate[1]);
        ((MCMoveStepTracker) sim.mcMoveTranslate[1].getTracker()).setNoisyAdjustment(true);
        ((MCMoveStepTracker) sim.mcMoveRotate[1].getTracker()).setNoisyAdjustment(true);

        MCMoveClusterConformationMDTest mcMove = new MCMoveClusterConformationMDTest(sim, space, temperature, f, l, true, true);
        MCMoveClusterConformationMDTest mcMove2 = new MCMoveClusterConformationMDTest(sim, space, temperature, f, l, true, false);
        mcMove2.setConformationsList(mcMove.getConformationsList());

        sim.integrators[0].getMoveManager().addMCMove(mcMove);
        sim.integrators[1].getMoveManager().addMCMove(mcMove2);

        // Here it defines the box length of every dimension
        double lb = sigmaTranslate * 3 + 2;
        sim.box[0].getBoundary().setBoxSize(Vector.of(new double[]{lb, lb, lb}));
        sim.box[1].getBoundary().setBoxSize(Vector.of(new double[]{lb, lb, lb}));

        if (false) {
            sim.box[0].getBoundary().setBoxSize(Vector.of(new double[]{lb, lb, lb}));
            sim.box[1].getBoundary().setBoxSize(Vector.of(new double[]{lb, lb, lb}));
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, "foo", 1);
            DisplayBox displayBox0 = simGraphic.getDisplayBox(sim.box[0]);
            DisplayBox displayBox1 = simGraphic.getDisplayBox(sim.box[1]);
            displayBox0.setShowBoundary(false);
            displayBox1.setShowBoundary(false);
            ((DisplayBoxCanvasG3DSys) displayBox0.canvas).setBackgroundColor(Color.WHITE);
            ((DisplayBoxCanvasG3DSys) displayBox1.canvas).setBackgroundColor(Color.WHITE);
            Box refBox = sim.getBox(0);
            Box targetBox = sim.getBox(1);

            simGraphic.getDisplayBox(refBox).setLabel("Reference-System Sampling");
            simGraphic.getDisplayBox(targetBox).setLabel("Target-System Sampling");
//            ((DisplayBoxCanvasG3DSys) displayBox0.canvas).setBackgroundColor(Color.WHITE);

            ColorScheme colorScheme0 = new ColorSchemeRandomByMolecule(sim, sim.box[0], sim.getRandom());
            ColorScheme colorScheme1 = new ColorSchemeRandomByMolecule(sim, sim.box[1], sim.getRandom());

            displayBox0.setColorScheme(colorScheme0);
            displayBox1.setColorScheme(colorScheme1);
            simGraphic.makeAndDisplayFrame();

            sim.setAccumulatorBlockSize(1000);
            sim.integratorOS.setNumSubSteps(1000);

            sim.getController().removeAction(sim.ai);
            sim.getController().addAction(new IAction() {
                public void actionPerformed() {
                    sim.initRefPref(null, 1000);
                    sim.equilibrate(null, 2000);

                    sim.ai.setMaxSteps(Long.MAX_VALUE);
                }
            });
            sim.getController().addAction(sim.ai);

            if ((Double.isNaN(sim.refPref) || Double.isInfinite(sim.refPref) || sim.refPref == 0)) {
                throw new RuntimeException("Oops");
            }
            return;
        }

        long t1 = System.currentTimeMillis();
        steps = steps / 1000;

        sim.initRefPref(null, steps / 20);
        sim.equilibrate(null, steps / 10);
//        sim.initRefPref(null, 5);
//        sim.equilibrate(null, 10);
        sim.ai.setMaxSteps(1000);
        sim.setAccumulatorBlockSize(steps);
        sim.integratorOS.setNumSubSteps(steps);


        sim.getController().actionPerformed();

        if ((Double.isNaN(sim.refPref) || Double.isInfinite(sim.refPref) || sim.refPref == 0)) {
            throw new RuntimeException("Oops");
        }

        long t2 = System.currentTimeMillis();

        System.out.println("final reference step frequency " + sim.integratorOS.getIdealRefStepFraction());
        System.out.println("actual reference step frequency " + sim.integratorOS.getRefStepFraction());

        sim.printResults(refIntegral);

        if (flex) {
            DataGroup allData = (DataGroup) sim.accumulators[1].getData();
            IData dataAvg = allData.getData(sim.accumulators[1].AVERAGE.index);
            IData dataErr = allData.getData(sim.accumulators[1].ERROR.index);
            IData dataCov = allData.getData(sim.accumulators[1].BLOCK_COVARIANCE.index);
            // we'll ignore block correlation -- whatever effects are here should be in the full target results
            int nTotal = (targetDiagrams.length + 2);
            double oVar = dataCov.getValue(nTotal * nTotal - 1);

            for (int i = 0; i < targetDiagrams.length; i++) {
                if (targetDiagramNumbers[i] < 0) {
                    System.out.print("diagram " + (-targetDiagramNumbers[i]) + ("bc "));
                } else {
                    System.out.print("diagram " + targetDiagramNumbers[i]);
//                    if (fTargetDiagramNumbers[i] != 0) {
//                        System.out.print(fTargetDiagramNumbers[i]);
//                    }

                    if (diagramFlexCorrection[i]) {
                        System.out.print("c");
                    }
                    System.out.print(" ");
                }
                // average is vi/|v| average, error is the uncertainty on that average
                // ocor is the correlation coefficient for the average and overlap values (vi/|v| and o/|v|)
                double ivar = dataCov.getValue((i + 1) * nTotal + (i + 1));
                double ocor = ivar * oVar == 0 ? 0 : dataCov.getValue(nTotal * (i + 1) + nTotal - 1) / Math.sqrt(ivar * oVar);
                System.out.print(String.format("average: %20.15e  error: %10.15e  ocor: %7.5f", dataAvg.getValue(i + 1), dataErr.getValue(i + 1), ocor));
                if (targetDiagrams.length > 1) {
                    System.out.print("  dcor:");
                    for (int j = 0; j < targetDiagrams.length; j++) {
                        if (i == j) continue;
                        double jvar = dataCov.getValue((j + 1) * nTotal + (j + 1));
                        double dcor = ivar * jvar == 0 ? 0 : dataCov.getValue((i + 1) * nTotal + (j + 1)) / Math.sqrt(ivar * jvar);
                        System.out.print(String.format(" %8.6f", dcor));
                    }
                }
                System.out.println();
            }
        }

        double AcceptProbTrans = sim.integrators[0].getMoveManager().getMCMoves().get(0).getTracker().acceptanceProbability();
        double AcceptProbRotate = sim.integrators[0].getMoveManager().getMCMoves().get(1).getTracker().acceptanceProbability();
        System.out.println();
        System.out.println("Reference System:");
        System.out.print(String.format("Acceptance Probability of McMoveTranslate: %f\n", AcceptProbTrans));
        System.out.print(String.format("Acceptance Probability of McMoveRotate: %f\n", AcceptProbRotate));

        double AcceptProbTrans2 = sim.mcMoveTranslate[1].getTracker().acceptanceProbability();
        double AcceptProbRotate2 = sim.mcMoveRotate[1].getTracker().acceptanceProbability();
        double AcceptProbConf2 = mcMove2.getTracker().acceptanceProbability();
        System.out.println("Target System:");
        System.out.print(String.format("Acceptance Probability of McMoveTranslate: %f\n", AcceptProbTrans2));
        System.out.print(String.format("Acceptance Probability of McMoveRotate: %f\n", AcceptProbRotate2));
        System.out.print(String.format("Acceptance Probability of McMoveConf: %f\n", AcceptProbConf2));
        System.out.println();
        System.out.println("time: " + (t2 - t1) / 1000.0);
    }

    /**
     * Inner class for parameters
     */
    public static class VirialHSParam extends ParameterBase {
        public int nPoints = 3;
        public long numSteps = 100000000;
        public static final int TREE = 0, CHAINS = 1;
        public int ref = CHAINS;
        public double chainFrac = 1;
        public int fv = 5;
        public int lv = 40;
    }

    public static String getSplitGraphString(Set<Graph> gSplit, VirialDiagrams flexDiagrams, boolean correction) {
        DeleteEdge edgeDeleter = new DeleteEdge();
        DeleteEdgeParameters ede = new DeleteEdgeParameters(flexDiagrams.eBond);
        boolean first = true;
        String str = "";
        for (Graph gs : gSplit) {
            if (VirialDiagrams.graphHasEdgeColor(gs, flexDiagrams.efbcBond)) {
                str += " " + gs.nodeCount() + "bc";
            } else {
                str += " " + gs.getStore().toNumberString();
                if (VirialDiagrams.graphHasEdgeColor(gs, flexDiagrams.eBond)) {
                    str += "p" + edgeDeleter.apply(gs, ede).getStore().toNumberString();
                }
            }
            if (first && correction) str += "c";
            first = false;
        }
        return str;
    }

}
