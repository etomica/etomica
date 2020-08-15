/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;

import etomica.action.activity.ActivityIntegrate2;
import etomica.box.Box;
import etomica.chem.elements.ElementSimple;
import etomica.data.histogram.HistogramSimple;
import etomica.graphics.ColorSchemeRandomByMolecule;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayBoxCanvasG3DSys;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorEvent;
import etomica.integrator.IntegratorListener;
import etomica.math.DoubleRange;
import etomica.math.SpecialFunctions;
import etomica.molecule.IMoleculeList;
import etomica.potential.IPotential;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphere;
import etomica.potential.Potential2SoftSpherical;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.util.random.RandomMersenneTwister;
import etomica.virial.*;
import etomica.virial.cluster.Standard;

import java.awt.*;
import java.io.File;
import java.util.Arrays;

/**
 * LJ simulation using Mayer sampling to evaluate cluster integrals.  Target
 * sampling weight is the umbrella sum of the derivative values.
 */
public class VirialLJDU {


    public static void main(String[] args) {

        VirialLJParam params = new VirialLJParam();
        if (args.length > 0) {
            params.writeRefPref = true;
        }
        ParseArgs.doParseArgs(params, args);
        if (args.length == 0) {
            params.nPoints = 5;
            params.nDer = 1;
            params.temperature = 1;
            params.numSteps = 10000000L;
            params.doHist = false;
            params.doChainRef = true;
            params.uWeights = null;
            params.writeRefPref = true;
            params.restartFilename = "foo";
        }

        runVirial(params);
    }

    public static void runVirial(VirialLJParam params) {

        if (params.sigmaHSRef == 0) {
            // below T=5, well is important and ref can sample a large HS
            if (params.temperature < 5) params.sigmaHSRef = 1.5;
                // above T=5, well is unimportant and LJ acts like SS.  range
                // continues to diminish as the temperature increases.
            else if (params.temperature < Double.POSITIVE_INFINITY) {
                params.sigmaHSRef = 0.925 * Math.pow(10 / params.temperature, 1.0 / 12.0);
                // sigma optimized for B7, but should be (approximately) universal
            } else {
                // SS with epsilon=4, temperature=1
                params.sigmaHSRef = 0.925 * Math.pow(10, 1.0 / 12.0);
            }
        }

        final int nPoints = params.nPoints;
        double temperature = params.temperature;
        if (temperature == Double.POSITIVE_INFINITY) temperature = 1;
        final int nDer = params.nDer;
        long steps = params.numSteps;
        double sigmaHSRef = params.sigmaHSRef;
        double refFrac = params.refFrac;
        boolean doHist = params.doHist;
        boolean doChainRef = params.doChainRef;
        long blockSize = params.blockSize;
        double[] uWeights = params.uWeights;
        boolean autoWeights = uWeights == null;
        if (autoWeights) {
            uWeights = new double[nDer + 1];
            uWeights[0] = 1;
        }

        double vhs = (4.0 / 3.0) * Math.PI * sigmaHSRef * sigmaHSRef * sigmaHSRef;
        final double HSBn = doChainRef ? SpecialFunctions.factorial(nPoints) / 2 * Math.pow(vhs, nPoints - 1) : Standard.BHS(nPoints, sigmaHSRef);
        System.out.println("sigmaHSRef: " + sigmaHSRef);
        System.out.println("B" + nPoints + "HS: " + HSBn);
        System.out.println("Lennard Jones overlap sampling B" + nPoints + " at T=" + temperature);
        Space space = Space3D.getInstance();

        MayerFunction fRefPos = new MayerFunction() {
            public void setBox(Box box) {
            }

            public IPotential getPotential() {
                return null;
            }

            public double f(IMoleculeList pair, double r2, double beta) {
                return r2 < sigmaHSRef * sigmaHSRef ? 1 : 0;
            }
        };

        MayerHardSphere fRef = new MayerHardSphere(sigmaHSRef);
        Potential2SoftSpherical pTarget = params.temperature < Double.POSITIVE_INFINITY ? new P2LennardJones(space) : new P2SoftSphere(space, 1, 4, 12);

        MayerGeneralSpherical fTarget = new MayerGeneralSpherical(pTarget);
        if (doChainRef) System.out.println("HS Chain reference");
        ClusterAbstract refCluster = doChainRef ? new ClusterChainHS(nPoints, fRefPos) : new ClusterWheatleyHS(nPoints, fRef);
        refCluster.setTemperature(temperature);

        final ClusterWheatleySoftDerivatives targetCluster = new ClusterWheatleySoftDerivatives(nPoints, fTarget, params.BDtol, nDer);
        targetCluster.setTemperature(temperature);

        if (blockSize == 0) blockSize = steps / 1000;
        System.out.println(steps + " steps (" + (steps / blockSize) + " blocks of " + blockSize + ")");

        ClusterMultivalueUmbrella targetUmbrella = new ClusterMultivalueUmbrella(targetCluster);
        targetUmbrella.setWeightCoefficients(uWeights);

        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space, new SpeciesSpheresMono(space, new ElementSimple("A")), nPoints, temperature, refCluster, targetCluster);
        if (params.randomSeeds == null) {
            System.out.println("random seeds: " + Arrays.toString(sim.getRandomSeeds()));
        } else {
            System.out.println("random seeds (set): " + Arrays.toString(params.randomSeeds));
            ((RandomMersenneTwister) sim.getRandom()).setSeedArray(params.randomSeeds);
        }
        sim.setSampleClusters(new ClusterWeight[]{new ClusterWeightAbs(refCluster), targetUmbrella});
        ClusterAbstract[] targetDiagrams = new ClusterAbstract[2 * (nDer + 1)];
        for (int m = 0; m <= nDer; m++) {
            ClusterWheatleySoftDerivatives.ClusterRetrievePrimes c = new ClusterWheatleySoftDerivatives.ClusterRetrievePrimes(targetCluster, m);
            targetDiagrams[m] = autoWeights ? new ClusterWeightAbs(c) : c;
        }
        for (int m = 0; m <= nDer; m++) {
            ClusterWheatleySoftDerivatives.ClusterRetrievePrimesBD c = new ClusterWheatleySoftDerivatives.ClusterRetrievePrimesBD(targetCluster, m);
            targetDiagrams[nDer + 1 + m] = autoWeights ? new ClusterWeightAbs(c) : c;
        }
        sim.setExtraTargetClusters(targetDiagrams);
        targetCluster.setBDAccFrac(params.BDAccFrac, sim.getRandom());
        sim.init();

        if (autoWeights) {
            String uWeightsFileName = params.writeRefPref ? "uWeights" + nPoints + "_" + temperature : null;
            double[] foo = sim.findUmbrellaWeights(uWeightsFileName, steps / 10);
            for (int i = 0; i < uWeights.length; i++) uWeights[i] = foo[i];

            targetUmbrella.setWeightCoefficients(uWeights);
            ((ClusterMultivalueUmbrella) sim.meters[0].getClusters()[1]).setWeightCoefficients(uWeights);
        } else {
            System.out.println("umbrella weights (set explicitly): " + Arrays.toString(uWeights));
        }

        long refBlocksize = blockSize;

        if (doChainRef) {
            sim.integrators[0].getMoveManager().removeMCMove(sim.mcMoveTranslate[0]);
            MCMoveClusterAtomHSChain mcMoveHSC = new MCMoveClusterAtomHSChain(sim.getRandom(), space, sigmaHSRef);
            sim.integrators[0].getMoveManager().addMCMove(mcMoveHSC);
            refBlocksize = Math.max(steps / (100 * 1000), 1);
            sim.accumulators[0].setBlockSize(refBlocksize);
        }

        int subSteps = 1000;
        sim.integratorOS.setNumSubSteps(subSteps);

        sim.integratorOS.setAggressiveAdjustStepFraction(true);

        if (false) {
            sim.box[0].getBoundary().setBoxSize(Vector.of(new double[]{10, 10, 10}));
            sim.box[1].getBoundary().setBoxSize(Vector.of(new double[]{10, 10, 10}));
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
            DisplayBox displayBox0 = simGraphic.getDisplayBox(sim.box[0]);
            DisplayBox displayBox1 = simGraphic.getDisplayBox(sim.box[1]);
//            displayBox0.setPixelUnit(new Pixel(300.0/size));
//            displayBox1.setPixelUnit(new Pixel(300.0/size));
            displayBox0.setShowBoundary(false);
            displayBox1.setShowBoundary(false);
            ((DisplayBoxCanvasG3DSys) displayBox0.canvas).setBackgroundColor(Color.WHITE);
            ((DisplayBoxCanvasG3DSys) displayBox1.canvas).setBackgroundColor(Color.WHITE);


            ColorSchemeRandomByMolecule colorScheme = new ColorSchemeRandomByMolecule(sim, sim.box[0], sim.getRandom());
            displayBox0.setColorScheme(colorScheme);
            colorScheme = new ColorSchemeRandomByMolecule(sim, sim.box[1], sim.getRandom());
            displayBox1.setColorScheme(colorScheme);
            simGraphic.makeAndDisplayFrame();

            sim.integratorOS.setNumSubSteps(1000);
            sim.setAccumulatorBlockSize(1000);

            // if running interactively, set filename to null so that it doens't read
            // (or write) to a refpref file
            sim.initRefPref(null, 10, false);
            sim.equilibrate(null, 20);
            sim.getController().addActivity(new ActivityIntegrate2(sim.integratorOS));
            if ((Double.isNaN(sim.refPref) || Double.isInfinite(sim.refPref) || sim.refPref == 0)) {
                throw new RuntimeException("Oops");
            }
            return;
        }

        long t1 = System.nanoTime();
        // if running interactively, don't use the file
        String refFileName = params.writeRefPref ? "refpref" + nPoints + "_" + temperature : null;
        // this will either read the refpref in from a file or run a short simulation to find it
        sim.initRefPref(refFileName, (steps / subSteps) / 20);

        // check for restart file and continue from there if possible
        File rf = params.restartFilename == null ? null : new File(params.restartFilename + ".sim");
        if (sim.getNumAlpha() == 1 && rf != null && rf.exists() && rf.length() > 0) {
            sim.setAccumulatorBlockSize(blockSize);
            sim.accumulators[0].setBlockSize(refBlocksize);
            long readSteps = sim.readRestart(params.restartFilename, false);
            steps -= readSteps;
        } else {
            // run another short simulation to find MC move step sizes and maybe narrow in more on the best ref pref
            // if it does continue looking for a pref, it will write the value to the file
            sim.equilibrate(refFileName, (steps / subSteps) / 10);
            System.out.println("equilibration finished");
            sim.setAccumulatorBlockSize(blockSize);
            if (doChainRef) sim.accumulators[0].setBlockSize(refBlocksize);
        }

        if (refFrac >= 0) {
            sim.integratorOS.setRefStepFraction(refFrac);
            sim.integratorOS.setAdjustStepFraction(false);
        }


        final HistogramSimple targHist = new HistogramSimple(200, new DoubleRange(-1, 4));
        IntegratorListener histListenerTarget = new IntegratorListener() {
            public void integratorStepStarted(IntegratorEvent e) {
            }

            public void integratorStepFinished(IntegratorEvent e) {
                CoordinatePairSet cPairs = sim.box[1].getCPairSet();
                for (int i = 0; i < nPoints; i++) {
                    for (int j = i + 1; j < nPoints; j++) {
                        double r2 = cPairs.getr2(i, j);
                        double r = Math.sqrt(r2);
                        if (r > 1) {
                            r = Math.log(r);
                        } else {
                            r -= 1;
                        }
                        targHist.addValue(r);
                    }
                }

            }

            public void integratorInitialized(IntegratorEvent e) {
            }
        };

        if (doHist) {
            System.out.println("collecting histograms");
            // only collect the histogram if we're forcing it to run the reference system
            sim.integrators[1].getEventManager().addListener(histListenerTarget);
        }


        IntegratorListener progressReport = new IntegratorListener() {

            public void integratorStepStarted(IntegratorEvent e) {
            }

            public void integratorStepFinished(IntegratorEvent e) {
                if (sim.integratorOS.getStepCount() % 100 != 0) return;
                System.out.print(sim.integratorOS.getStepCount() + " steps: ");
                double[] ratioAndError = sim.dvo.getAverageAndError();
                System.out.println("abs average: " + ratioAndError[0] * HSBn + ", error: " + ratioAndError[1] * HSBn);
            }

            public void integratorInitialized(IntegratorEvent e) {
            }

        };
        if (false) {
            sim.integratorOS.getEventManager().addListener(progressReport);
        }

        if (params.restartFilename != null) {
            sim.setWriteRestart(params.restartFilename, false);
        }

        sim.integratorOS.setNumSubSteps((int) blockSize);
        for (int i = 0; i < 2; i++) {
            if (i > 0 || !doChainRef) System.out.println("MC Move step sizes " + sim.mcMoveTranslate[i].getStepSize());
        }
sim.getController().runActivityBlocking(new ActivityIntegrate2(sim.integratorOS), steps / blockSize);
        long t2 = System.nanoTime();

        if (doHist) {
            double[] xValues = targHist.xValues();
            double[] h = targHist.getHistogram();
            for (int i = 0; i < xValues.length; i++) {
                if (!Double.isNaN(h[i])) {
                    double r = xValues[i];
                    double y = h[i];
                    if (r < 0) r += 1;
                    else {
                        r = Math.exp(r);
                        y /= r;
                    }
                    System.out.println(r + " " + y);
                }
            }
        }

        System.out.println("final reference step frequency " + sim.integratorOS.getIdealRefStepFraction());
        System.out.println("actual reference step frequency " + sim.integratorOS.getRefStepFraction());

        String[] extraNames = new String[2 * (nDer + 1)];
        for (int i = 0; i <= nDer; i++) {
            extraNames[i] = "derivative " + i;
        }
        for (int i = 0; i <= nDer; i++) {
            extraNames[nDer + 1 + i] = "derivativeBD " + i;
        }
        sim.printResults(HSBn, extraNames);

        double timeBD = targetCluster.getTimeBD();
        System.out.println("timeBDfrac: " + timeBD / (t2 - t1) * 1e9);
        System.out.println("time: " + (t2 - t1) / 1e9);
        if (params.restartFilename != null) {
            new File(params.restartFilename + ".sim").delete();
            new File(params.restartFilename + ".ref").delete();
            new File(params.restartFilename + ".target").delete();
        }
    }

    /**
     * Inner class for parameters
     */
    public static class VirialLJParam extends ParameterBase {
        public int nPoints = 4;
        public double temperature = 1;
        public int nDer = 2;
        public long numSteps = 10000000;
        public double sigmaHSRef = 0;
        public boolean writeRefPref = false;
        public double refFrac = -1;
        public boolean doHist = false;
        public boolean doChainRef = true;
        public long blockSize = 0;
        public double BDAccFrac = 0.1;
        public double BDtol = 1e-12;
        public double[] uWeights = null;
        public int[] randomSeeds = null;
        public String restartFilename = null;
    }
}
