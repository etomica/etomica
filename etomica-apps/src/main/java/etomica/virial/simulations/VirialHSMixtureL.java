/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.DiameterHash;
import etomica.atom.IAtom;
import etomica.chem.elements.ElementSimple;
import etomica.graphics.*;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.species.SpeciesSpheresMono;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.virial.ClusterWheatleyHS;
import etomica.virial.MCMoveClusterAtomMulti;
import etomica.virial.MayerFunction;
import etomica.virial.MayerHSMixture;

import javax.swing.*;
import java.awt.*;

/**
 * Calculation for virial coefficients of hard sphere mixture
 *
 * @author Andrew
 */
public class VirialHSMixtureL {

    public static void main(String[] args) {

        VirialHSParam params = new VirialHSParam();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        } else {
            params.nPoints = 6;
            params.numSteps = 1000000000L;
            params.q = 0.2;
            params.nonAdd = 0.0;
            params.cff = 4;
            params.D = 3;
            params.targetL = 2;
        }
        final int nSmall = params.cff;
        final int nPoints = params.nPoints;
        long steps = params.numSteps;
        final double q = params.q;
        System.out.println("   q: " + q);
        final double sigmaB = 1;
        final double sigmaS = q * sigmaB;
        final double[][] nonAdd = new double[2][2];
        nonAdd[0][1] = nonAdd[1][0] = params.nonAdd;
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

        double refL = nBig + nSmall * q;
        BoundaryRectangularPeriodic boundaryRef = new BoundaryRectangularPeriodic(space, refL);
        BoundaryRectangularPeriodic boundaryTarget = new BoundaryRectangularPeriodic(space, params.targetL);
        MayerHSMixture fTarget = MayerHSMixture.makeTargetF(space, nPoints, pairSigma, boundaryTarget);
        MayerFunction fRef = MayerHSMixture.makeTargetF(space, nPoints, pairSigma, boundaryRef);

        final ClusterWheatleyHS targetCluster = new ClusterWheatleyHS(nPoints, fTarget);
        targetCluster.setTemperature(1.0);
        final ClusterWheatleyHS refCluster = new ClusterWheatleyHS(nPoints, fRef);
        refCluster.setTemperature(1.0);
        final double refIntegral = Math.pow(params.targetL / refL, params.D * (nPoints - 1));
        System.out.println("reference integral: " + refIntegral);
        refCluster.setTemperature(1.0);

        long blockSize = steps / 1000;
        System.out.println(steps + " steps (" + (steps / blockSize) + " blocks of " + blockSize + ")");
        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space, new SpeciesSpheresMono(space, new ElementSimple("A")), nPoints, 1, refCluster, targetCluster);
        sim.setBoxLengths(refL, params.targetL);
        sim.init();
        sim.integratorOS.setDoAdjustOnTime(true);

        ((MCMoveClusterAtomMulti) sim.mcMoveTranslate[0]).setStartAtom(0);
        ((MCMoveClusterAtomMulti) sim.mcMoveTranslate[0]).setDoImposePBC(true);
        ((MCMoveClusterAtomMulti) sim.mcMoveTranslate[1]).setStartAtom(0);
        ((MCMoveClusterAtomMulti) sim.mcMoveTranslate[1]).setDoImposePBC(true);

        int subSteps = 1000;
        sim.integratorOS.setNumSubSteps(subSteps);

        sim.integratorOS.setAggressiveAdjustStepFraction(true);

        if (true) {
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
            DisplayBox displayBox0 = simGraphic.getDisplayBox(sim.box[0]);
            DisplayBox displayBox1 = simGraphic.getDisplayBox(sim.box[1]);

            DiameterHash dh = new DiameterHash() {
                @Override
                public double getDiameter(IAtom atom) {
                    int idx = atom.getLeafIndex();
                    if (idx < nBig) return 1;
                    return q;
                }
            };
            displayBox0.setDiameterHash(dh);
            displayBox1.setDiameterHash(dh);


            ColorScheme colorScheme = new ColorScheme() {

                public Color getAtomColor(IAtom a) {
                    float b = a.getLeafIndex() / ((float) nPoints);
                    float r = 1.0f - b;
                    int idx = a.getLeafIndex();
                    return new Color(r, 0f, b, idx < nBig ? 0.1f : 1.0f);
                }
            };
            displayBox0.setColorScheme(colorScheme);
            displayBox1.setColorScheme(colorScheme);
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
        for (int i = 0; i < 2; i++) {
            System.out.println("MC Move step size (" + i + ") " + sim.mcMoveTranslate[i].getStepSize());
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
        public double targetL = Double.NaN;
        public double q = 0.1;
        public int cff = 0;
        public double nonAdd = 0;
        public int D = 3;
        public double refFrac = -1;
        public boolean writeRefPref = false;
    }

}
