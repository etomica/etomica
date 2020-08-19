/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.IAtom;
import etomica.chem.elements.ElementSimple;
import etomica.graphics.*;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.SpeciesSpheresMono;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.virial.*;

import javax.swing.*;
import java.awt.*;

/**
 * Calculation for virial coefficients of hard sphere mixture
 *
 * @author Andrew
 */
public class VirialHSMixtureOverlap {

    public static void main(String[] args) {

        VirialHSParam params = new VirialHSParam();
        if (args.length > 0) {
            params.writeRefPref = true;
            ParseArgs.doParseArgs(params, args);
        } else {
            params.nPoints = 8;
            params.numSteps = 100000000L;
            params.ref = RefChoice.CHAIN_TREE;
            params.chainFrac = 0.5;
            params.q = 0.2;
            params.nonAdd = 0.0;
            params.cff = 6;
            params.D = 3;
        }
        final int nSmall = params.cff;
        final int nPoints = params.nPoints;
        long steps = params.numSteps;
        final double q = params.q;
        final RefChoice ref = params.ref;
        System.out.println("   q: " + q);
        final double sigmaB = 1;
        final double sigmaS = q * sigmaB;
        final double[][] nonAdd = new double[2][2];
        nonAdd[0][1] = nonAdd[1][0] = params.nonAdd;
        final double chainFrac = params.chainFrac;
        final int D = params.D;
        final double refFrac = params.refFrac;
        System.out.println("Number of points: " + nPoints);
        System.out.println("Dimensions: " + D);
        System.out.println("nonAdd: " + params.nonAdd);
        Space space = Space.getInstance(D);

        int nBig = nPoints - nSmall;
        System.out.println("B" + nBig + nSmall);
        int[] nTypes = new int[]{nBig, nSmall};
        double[] sigma = new double[]{sigmaB, sigmaS};
        double[][] pairSigma = new double[nPoints][nPoints];

        int iii = 0;
        for (int i = 0; i < nTypes.length; i++) {
            for (int ii = 0; ii < nTypes[i]; ii++) {
                int jjj = 0;
                for (int j = 0; j <= i; j++) {
                    for (int jj = 0; jj <= ((i == j) ? ii : (nTypes[j] - 1)); jj++) {
                        double x = nonAdd[i][j];
                        pairSigma[jjj][iii] = pairSigma[iii][jjj] = (1 + x) * (sigma[i] + sigma[j]) / 2;
                        jjj++;
                    }
                }
                iii++;
            }
        }

        MayerHSMixture fTarget = MayerHSMixture.makeTargetF(space, nPoints, pairSigma, null);
        MayerFunction fRefPos = MayerHSMixture.makeReferenceF(space, nPoints, pairSigma, null);

        final ClusterAbstract targetCluster = new ClusterWheatleyHS(nPoints, fTarget);
        targetCluster.setTemperature(1.0);
        final ClusterWeightUmbrella refCluster;
        double ri;
        if (ref == RefChoice.TREE) {
            System.out.println("using a tree reference");
            ClusterSinglyConnected ct = new ClusterSinglyConnected(nPoints, fRefPos);
            refCluster = new ClusterWeightUmbrella(new ClusterAbstract[]{ct});
            refCluster.setWeightCoefficients(new double[]{1.0 / ct.numDiagrams()});
            ri = 1;
        } else if (ref == RefChoice.CHAINS) {
            System.out.println("using a chain reference");
            ClusterChainHS cc = new ClusterChainHS(nPoints, fRefPos);
            refCluster = new ClusterWeightUmbrella(new ClusterAbstract[]{cc});
            refCluster.setWeightCoefficients(new double[]{1.0 / cc.numDiagrams()});
            ri = 1;
        } else if (ref == RefChoice.CHAIN_TREE) {
            System.out.println("using a chain/tree reference (" + chainFrac + " chains)");
            ClusterChainHS cc = new ClusterChainHS(nPoints, fRefPos);
            ClusterSinglyConnected ct = new ClusterSinglyConnected(nPoints, fRefPos);
            refCluster = new ClusterWeightUmbrella(new ClusterAbstract[]{cc, ct});
            refCluster.setWeightCoefficients(new double[]{chainFrac / cc.numDiagrams(), (1 - chainFrac) / ct.numDiagrams()});
            ri = 1;
        } else {
            throw new RuntimeException();
        }

        final double refIntegral = ri;
        System.out.println("reference integral: " + refIntegral);
        refCluster.setTemperature(1.0);

        long blockSize = steps / 1000;
        System.out.println(steps + " steps (" + (steps / blockSize) + " blocks of " + blockSize + ")");
        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space, new SpeciesSpheresMono(space, new ElementSimple("A")), nPoints, 1, refCluster, targetCluster);
        sim.init();
        sim.integratorOS.setDoAdjustOnTime(true);

        sim.integrators[0].getMoveManager().removeMCMove(sim.mcMoveTranslate[0]);
        if (ref == RefChoice.TREE) {
            MCMoveClusterAtomHSTreeMix mcMoveHS = new MCMoveClusterAtomHSTreeMix(sim.getRandom(), space, pairSigma);
            sim.integrators[0].getMoveManager().addMCMove(mcMoveHS);
        } else if (ref == RefChoice.CHAINS) {
            MCMoveClusterAtomHSChainMix mcMoveHS = new MCMoveClusterAtomHSChainMix(sim.getRandom(), space, pairSigma);
            sim.integrators[0].getMoveManager().addMCMove(mcMoveHS);
        } else if (ref == RefChoice.CHAIN_TREE) {
            MCMoveClusterAtomHSTreeMix mcMoveHST = new MCMoveClusterAtomHSTreeMix(sim.getRandom(), space, pairSigma);
            sim.integrators[0].getMoveManager().addMCMove(mcMoveHST);
            sim.integrators[0].getMoveManager().setFrequency(mcMoveHST, 1 - chainFrac);
            MCMoveClusterAtomHSChainMix mcMoveHSC = new MCMoveClusterAtomHSChainMix(sim.getRandom(), space, pairSigma);
            sim.integrators[0].getMoveManager().addMCMove(mcMoveHSC);
            sim.integrators[0].getMoveManager().setFrequency(mcMoveHSC, chainFrac);
        }

        int subSteps = 1000;
        sim.integratorOS.setNumSubSteps(subSteps);

        sim.integratorOS.setAggressiveAdjustStepFraction(true);

        if (false) {
            sim.box[0].getBoundary().setBoxSize(Vector.of(new double[]{10, 10, 10}));
            sim.box[1].getBoundary().setBoxSize(Vector.of(new double[]{10, 10, 10}));
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
            DisplayBox displayBox0 = simGraphic.getDisplayBox(sim.box[0]);
            displayBox0.setShowBoundary(false);
            ((DisplayBoxCanvasG3DSys) displayBox0.canvas).setBackgroundColor(Color.WHITE);
//            displayBox0.setPixelUnit(new Pixel(300.0/size));
            DisplayBox displayBox1 = simGraphic.getDisplayBox(sim.box[1]);
            displayBox1.setShowBoundary(false);
            ((DisplayBoxCanvasG3DSys) displayBox1.canvas).setBackgroundColor(Color.WHITE);


            ColorScheme colorScheme = new ColorScheme() {

                public Color getAtomColor(IAtom a) {
                    float b = a.getLeafIndex() / ((float) nPoints);
                    float r = 1.0f - b;
                    return new Color(r, 0f, b);
                }
            };
            displayBox0.setColorScheme(colorScheme);
            simGraphic.makeAndDisplayFrame();

            sim.setAccumulatorBlockSize(1000);

            final JPanel panelParentGroup = new JPanel(new GridBagLayout());
            GridBagConstraints gbc = new GridBagConstraints();
            gbc.gridx = 0;
            gbc.gridwidth = 2;
            final DisplayTextBox stepsBox = new DisplayTextBox();
            stepsBox.setLabel("steps");
            panelParentGroup.add(stepsBox.graphic(), gbc);

            final DisplayTextBox averageBox = new DisplayTextBox();
            averageBox.setLabel("Average");
            final DisplayTextBox errorBox = new DisplayTextBox();
            errorBox.setLabel("Error");
            JLabel jLabelPanelParentGroup = new JLabel("B" + nPoints);
            gbc.gridy = 1;
            panelParentGroup.add(jLabelPanelParentGroup, gbc);
            gbc.gridy = 2;
            gbc.gridwidth = 1;
            panelParentGroup.add(averageBox.graphic(), gbc);
            gbc.gridx = 1;
            panelParentGroup.add(errorBox.graphic(), gbc);
            simGraphic.getPanel().controlPanel.add(panelParentGroup, SimulationPanel.getVertGBC());

            return;
        }

        long t1 = System.currentTimeMillis();
        // if running interactively, don't use the file
        String refFileName = params.writeRefPref ? "refpref" + nPoints : null;
        // this will either read the refpref in from a file or run a short simulation to find it
        sim.initRefPref(refFileName, (steps / subSteps) / 20);
        // run another short simulation to find MC move step sizes and maybe narrow in more on the best ref pref
        // if it does continue looking for a pref, it will write the value to the file
        sim.equilibrate(refFileName, (steps / subSteps) / 10);
ActivityIntegrate ai = new ActivityIntegrate(sim.integratorOS, steps / blockSize);
System.out.println("equilibration finished");

        if (refFrac >= 0) {
            sim.integratorOS.setRefStepFraction(refFrac);
            sim.integratorOS.setAdjustStepFraction(false);
        }

        sim.integratorOS.setNumSubSteps((int) blockSize);
        sim.setAccumulatorBlockSize(blockSize);
        sim.accumulators[0].setBlockSize(1);
        for (int i = 1; i < 2; i++) {
            System.out.println("MC Move step sizes " + sim.mcMoveTranslate[i].getStepSize());
        }
sim.getController().runActivityBlocking(ai);
        long t2 = System.currentTimeMillis();

        System.out.println("final reference step fraction " + sim.integratorOS.getIdealRefStepFraction());
        System.out.println("actual reference step fraction " + sim.integratorOS.getRefStepFraction());
        System.out.println("time reference step fraction " + sim.integratorOS.getRefTimeFraction());

        sim.printResults(refIntegral, null);

        System.out.println("time:" + (t2 - t1) / 1000.0);
    }

    enum RefChoice {
        TREE, CHAINS, CHAIN_TREE
    }

    /**
     * Inner class for parameters
     */
    public static class VirialHSParam extends ParameterBase {
        public int nPoints = 3;
        public long numSteps = 100000000;
        public RefChoice ref = RefChoice.CHAIN_TREE;
        public double chainFrac = 0.5;
        public double q = 0.1;
        public int cff = 0;
        public double nonAdd = 0;
        public int D = 3;
        public double refFrac = -1;
        public boolean writeRefPref = false;
    }

}
