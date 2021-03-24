/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;

import etomica.atom.DiameterHash;
import etomica.atom.IAtom;
import etomica.chem.elements.ElementSimple;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.IData;
import etomica.data.types.DataGroup;
import etomica.graphics.*;
import etomica.integrator.mcmove.MCMoveStepTracker;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.species.ISpecies;
import etomica.species.SpeciesSpheresMono;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.virial.*;

import javax.swing.*;
import java.awt.*;

/**
 * Calculation for virial coefficients of hard sphere mixture
 *
 * @author Arpit
 */

public class VirialHSMixturePercolation {
    public static void main(String[] args) {

        VirialHSParam params = new VirialHSParam();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        } else {
            params.nPoints = 6;
            params.numSteps = 10000000L;
            params.q = 0.2;
            params.nonAdd = 0.0;
            params.cff = 4;
            params.D = 3;
            params.L = 1.8;
            params.order = "BSSBSS";
        }
        final int nSmall = params.cff;
        final int nPoints = params.nPoints;
        long steps = params.numSteps;
        final double q = params.q;
//        System.out.println("   q: " + q);
        final double sigmaB = 1;
        final double sigmaS = q * sigmaB;
        final double[][] nonAdd = new double[2][2];
        nonAdd[0][1] = nonAdd[1][0] = params.nonAdd;
        final int D = params.D;
        final double L = params.L;
//        System.out.println("Number of points: " + nPoints);
//        System.out.println("Dimensions: " + D);
//        System.out.println("nonAdd: " + params.nonAdd);
        Space space = Space.getInstance(D);
        long t1 = System.currentTimeMillis();

        int nBig = nPoints - nSmall;
//        System.out.println("B" + nBig + nSmall);
        int[] nTypes = new int[]{nBig, nSmall};
//        double[] sigma = new double[]{sigmaB, sigmaS};
//        double[] sigma = new double[]{sigmaB, sigmaS, sigmaS, sigmaB, sigmaS, sigmaS};
        double[] sigma = new double[nPoints];
        for (int i = 0; i < nPoints; i++) {
            if (params.order.charAt(i) == 'B')
                sigma[i]=sigmaB;
            else
                sigma[i]=sigmaS;
//            System.out.print(sigma[i]+" ");
        }
        double[][] pairSigma = new double[nPoints][nPoints];

//        int iii = 0;
//        for (int i = 0; i < nTypes.length; i++) {
//            for (int ii = 0; ii < nTypes[i]; ii++) {
//                int jjj = 0;
//                for (int j = 0; j <= i; j++) {
//                    for (int jj = 0; jj <= ((i == j) ? ii : (nTypes[j] - 1)); jj++) {
//                        double x = nonAdd[i][j];
//                        pairSigma[jjj][iii] = pairSigma[iii][jjj] = (1 + x) * (sigma [i] + sigma[j]) / 2;
//                        jjj++;
//                    }
//                }
//                iii++;
//            }
//        }
        for (int i = 0; i < sigma.length; i++) {
//            System.out.print(sigma[i]+" ");
            for (int j = 0; j < i+1; j++) {
                pairSigma[i][j] = pairSigma[j][i] = (sigma[i] + sigma[j]) / 2;
            }
        }
//        System.out.println();

        BoundaryRectangularPeriodic boundaryRef = new BoundaryRectangularPeriodic(space, params.L);
        MayerHSMixture fTarget = MayerHSMixture.makeTargetF(space, nPoints, pairSigma, boundaryRef);
        final ClusterAbstract refCluster = new ClusterConstant(nPoints, 1);
        final ClusterWheatleyHS subCluster = new ClusterWheatleyHS(nPoints, fTarget);
        subCluster.setTemperature(1.0);
        final ClusterChainSpanBox targetCluster = new ClusterChainSpanBox(nPoints, pairSigma, subCluster);
        targetCluster.setTemperature(1.0);
        final double refIntegral = 1.0;
//        System.out.println("reference integral: " + refIntegral);

        ClusterAbstract[] targetDiagrams = new ClusterAbstract[]{targetCluster};

        long blockSize = steps / 1000;
//        System.out.println(steps + " steps (" + (steps / blockSize) + " blocks of " + blockSize + ")");

        final SimulationVirial sim = new SimulationVirial(space, new ISpecies[]{new SpeciesSpheresMono(space, new ElementSimple("A"))}, new int[]{nPoints}, 1.0, ClusterWeightAbs.makeWeightCluster(targetCluster), refCluster, targetDiagrams);
//        int[] seeds = new int[]{1};
//        sim.setSeeds(seeds);
        if (L > 0 && L < Double.POSITIVE_INFINITY) sim.setBoxLength(L);
        sim.init();
        double factor = sim.boxLength/(sigmaB*nBig + sigmaS*nSmall);
        sim.box.getLeafList().get(0).getPosition().setX(0, -L/2 + factor*0.5*pairSigma[0][0]);
        for (int i=1; i<sim.box.getLeafList().size(); i++){
            sim.box.getLeafList().get(i).getPosition().setX(0, sim.box.getLeafList().get(i-1).getPosition().getX(0) + factor*pairSigma[i][i-1]);
        }
        MeterPercolation meter = new MeterPercolation(sim.box);
        sim.setMeter(meter);
        AccumulatorAverageFixed accumulator = new AccumulatorAverageFixed(params.numSteps/100);
        sim.setAccumulator(accumulator);
        accumulator.setPushInterval(1000000);

        sim.integrator.getMoveManager().removeMCMove(sim.mcMoveTranslate);
        MCMoveClusterFixCOM mcMoveClusterFixCOM = new MCMoveClusterFixCOM(sim.getRandom(), space);
        sim.integrator.getMoveManager().addMCMove(mcMoveClusterFixCOM);
        sim.integrator.getMoveManager().setEquilibrating(true);

        if(false){
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
            DisplayBox displayBox0 = simGraphic.getDisplayBox(sim.box);
            displayBox0.setDiameterHash(new DiameterHash() {
                @Override
                public double getDiameter(IAtom atom) {
                    return atom.getLeafIndex() < nBig ? 1 : q;
                }
            });
//            displayBox0.setPixelUnit(new Pixel(300.0/size));
            if (L == Double.POSITIVE_INFINITY || L == 0) {
                displayBox0.setShowBoundary(false);
                ((DisplayBoxCanvasG3DSys) displayBox0.canvas).setBackgroundColor(Color.WHITE);
            }

            ColorScheme colorScheme = new ColorScheme() {

                public Color getAtomColor(IAtom a) {
                    float b = a.getLeafIndex() / ((float) nPoints);
                    float r = 1.0f - b;
                    return new Color(r, 0f, b);
                }
            };
            displayBox0.setColorScheme(colorScheme);
            simGraphic.makeAndDisplayFrame();
            return;
        }

        sim.ai.setMaxSteps(steps/10);
        sim.getController().actionPerformed();
        sim.getController().reset();
        accumulator.reset();

        sim.ai.setMaxSteps(steps);
        sim.getController().actionPerformed();
//        System.out.print(mcMoveClusterFixCOM.getStepSize());

        DataGroup allYourBase = (DataGroup) accumulator.getData();
        IData averageData = allYourBase.getData(accumulator.AVERAGE.index);
        IData errorData = allYourBase.getData(accumulator.ERROR.index);

//        System.out.println();

        for (int i=0; i<sim.box.getLeafList().size(); i++){
//            System.out.print(i);
            for (int j=0; j<sim.getSpace().D(); j++){
                double avg = Math.sqrt(averageData.getValue(i*sim.getSpace().D() + j));
                double err = (errorData.getValue(i*sim.getSpace().D() + j))/(2*avg);
                System.out.print(String.format("  %20.15e %9.4e ", avg, err));
//                System.out.print(String.format("%20.15e ", avg));
            }
            System.out.println();
        }

        long t2 = System.currentTimeMillis();
        System.out.println("time:" + (t2 - t1) / 1000.0);
//        System.out.println("#################################");

    }
    /**
     * Inner class for parameters
     */
    public static class VirialHSParam extends ParameterBase {
        public int nPoints = 3;
        public long numSteps = 100000000;
        public double L = Double.NaN;
        public double q = 0.1;
        public int cff = 0;
        public double nonAdd = 0;
        public int D = 3;
        public String order = "BBSSSS";
    }

}

