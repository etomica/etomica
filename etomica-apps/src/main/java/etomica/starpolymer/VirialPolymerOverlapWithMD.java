/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.starpolymer;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.chem.elements.ElementSimple;
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
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.virial.*;
import etomica.virial.cluster.VirialDiagrams;
import etomica.virial.simulations.SimulationVirialOverlap2;

import java.awt.*;
import java.util.Set;

/**
 * Calculation for virial coefficients of hard spheres
 * Edited for calculating protein 1HRB
 * <p>
 * CHANGES: change for star-polymers. May 7, 2018
 */
public class VirialPolymerOverlapWithMD {

    public static void main(String[] args) {

        VirialHSParam params = new VirialHSParam();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        } else {
            params.nPoints = 2;
            params.numSteps = 100000L;
            params.ref = VirialHSParam.CHAINS;
            params.fv = 3;
            params.lv = 34;
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
        System.out.println(steps + " steps");

        ClusterAbstract targetCluster = new ClusterWheatleyHS(nPoints, fRef);
        targetCluster.setTemperature(temperature);

        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space, species, nPoints, temperature, refCluster, targetCluster);
//        sim.setRandom(new RandomMersenneTwister(2));
        sim.init();

        if (ref == VirialHSParam.CHAINS) {
            sim.integrators[0].getMoveManager().removeMCMove(sim.mcMoveTranslate[0]);
            sim.integrators[0].getMoveManager().addMCMove(new MCMoveClusterMoleculeHSChain(sim.getRandom(), space, refDiameter));
            sim.accumulators[0].setBlockSize(1);
        }

        ((MCMoveClusterRotateMoleculeMulti) sim.mcMoveRotate[0]).setPositionDefinition(new MoleculePositionCOM(space));
//        sim.integrators[1].getMoveManager().removeMCMove(sim.mcMoveRotate[1]);
        ((MCMoveStepTracker) sim.mcMoveTranslate[1].getTracker()).setNoisyAdjustment(true);
        ((MCMoveStepTracker) sim.mcMoveRotate[1].getTracker()).setNoisyAdjustment(true);
//        sim.integrators[1].getMoveManager().setFrequency();

        MCMoveClusterConformationMDTest mcMove = new MCMoveClusterConformationMDTest(sim, space, temperature, f, l, true, true);
        System.exit(0);
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

            sim.initRefPref(null, 1000, false);
            sim.equilibrate(null, 2000);
            sim.getController().addActivity(new ActivityIntegrate(sim.integratorOS));

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
        sim.setAccumulatorBlockSize(steps);
        sim.integratorOS.setNumSubSteps(steps);
sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integratorOS), 1000);

        if ((Double.isNaN(sim.refPref) || Double.isInfinite(sim.refPref) || sim.refPref == 0)) {
            throw new RuntimeException("Oops");
        }

        long t2 = System.currentTimeMillis();

        System.out.println("final reference step frequency " + sim.integratorOS.getIdealRefStepFraction());
        System.out.println("actual reference step frequency " + sim.integratorOS.getRefStepFraction());

        sim.printResults(refIntegral);

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
