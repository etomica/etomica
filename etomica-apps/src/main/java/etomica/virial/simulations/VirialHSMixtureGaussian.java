package etomica.virial.simulations;

import etomica.atom.DiameterHash;
import etomica.atom.IAtom;
import etomica.chem.elements.ElementSimple;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.DataPumpListener;
import etomica.data.IData;
import etomica.data.types.DataGroup;
import etomica.graphics.*;
import etomica.integrator.mcmove.MCMoveStepTracker;
import etomica.math.numerical.ArrayReader1D;
import etomica.normalmode.ArrayReader2D;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.species.SpeciesSpheresMono;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.util.random.RandomMersenneTwister;
import etomica.virial.*;

import javax.swing.*;
import java.awt.*;
import java.io.FileWriter;
import java.io.IOException;

public class VirialHSMixtureGaussian {
    public static void main(String[] args) throws IOException {

        VirialHSMixtureGaussian.VirialHSParam params = new VirialHSMixtureGaussian.VirialHSParam();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        } else {
            params.nPoints = 6;
            params.numSteps = 10000000L;
            params.q = 0.2;
            params.nonAdd = 0.0;
            params.cff = 4;
            params.D = 3;
            params.targetL = 1.8;
            params.targetAcceptance = 0.2;
            params.order = "BSSBSS";
//            params.order = new char[]{'B', 'B', 'S', 'S', 'S', 'S'};
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
        final double targetAcceptance = params.targetAcceptance;
        System.out.println("Number of points: " + nPoints);
        System.out.println("Dimensions: " + D);
        System.out.println("nonAdd: " + params.nonAdd);
        System.out.println("targetAcceptance: " + targetAcceptance);
        Space space = Space.getInstance(D);

        int nBig = nPoints - nSmall;
        System.out.println("B" + nBig + nSmall);
        int[] nTypes = new int[]{nBig, nSmall};
//        double[] sigma = new double[]{sigmaB, sigmaS};
//        System.out.println("1"+params.order);
        double[] sigma = new double[nPoints];
        for (int i = 0; i < nPoints; i++) {
            if (params.order.charAt(i) == 'B')
                sigma[i]=sigmaB;
            else
                sigma[i]=sigmaS;
//            System.out.print(sigma[i]+" ");
        }
//        System.out.println("2"+params.order);
        double[][] pairSigma = new double[nPoints][nPoints];
//        double[][] standardDeviation = new double[nPoints][D];
        double [][] meanPosition = new double[nPoints][D];
//        FileWriter writer = new FileWriter("test.txt");
//        for (int i = 0; i < nPoints; i++) {
//            writer.write(i+"\n");
//        }
//        writer.close();
        double[][] standardDeviation = ArrayReader1D.getFromFile("StandardDeviationB"+nBig+nSmall+"L"+params.targetL+params.order+".txt");
        for (int i = 0; i < nPoints; ++i) {
            for (int j = 0; j < D; ++j) {
                System.out.print(standardDeviation[i][j]+" ");
            }
            System.out.println();
        }
//        standardDeviation[0][0] = 2.665100684371301e-03;
//        standardDeviation[0][1] = 1.188494201406047e-02;
//        standardDeviation[0][2] = 1.188410178860041e-02;
//        standardDeviation[1][0] = 2.663992573559614e-03;
//        standardDeviation[1][1] = 1.189003484762412e-02;
//        standardDeviation[1][2] = 1.189065780008437e-02;
//        standardDeviation[2][0] = 2.521717082184365e-03;
//        standardDeviation[2][1] = 5.605474215964236e-03;
//        standardDeviation[2][2] = 5.600831778281200e-03;
//        standardDeviation[3][0] = 2.520623180812475e-03;
//        standardDeviation[3][1] = 5.610869606540283e-03;
//        standardDeviation[3][2] = 5.610071243253531e-03;

//        int iii = 0;
//        for (int i = 0; i < nTypes.length; i++) {
//            for (int ii = 0; ii < nTypes[i]; ii++) {
//                int jjj = 0;
//                for (int j = 0; j <= i; j++) {
//                    for (int jj = 0; jj <= ((i == j) ? ii : (nTypes[j] - 1)); jj++) {
//                        double x = nonAdd[i][j];
//                        pairSigma[jjj][iii] = pairSigma[iii][jjj] = (1 + x) * (sigma[i] + sigma[j]) / 2;
//                        jjj++;
//                    }
//                }
//                iii++;
//            }
//        }
        for (int i = 0; i < sigma.length; i++) {
            for (int j = 0; j < i+1; j++) {
                pairSigma[i][j] = pairSigma[j][i] = (sigma[i] + sigma[j]) / 2;
//                System.out.println("pairSigma["+i+"]["+j+"]"+pairSigma[i][j]+" pairSigma["+j+"]["+i+"]"+pairSigma[j][i]);
            }
        }

//        double refL = nBig + nSmall * q;
        BoundaryRectangularPeriodic boundaryRef = new BoundaryRectangularPeriodic(space, params.targetL);
        MayerHSMixture fTarget = MayerHSMixture.makeTargetF(space, nPoints, pairSigma, boundaryRef);
        MayerFunction fRef = MayerHSMixture.makeTargetF(space, nPoints, pairSigma, boundaryRef);

        final ClusterWheatleyHS subCluster = new ClusterWheatleyHS(nPoints, fTarget);
        subCluster.setTemperature(1.0);
        final ClusterChainSpanBox targetCluster = new ClusterChainSpanBox(nPoints, pairSigma, subCluster);
        targetCluster.setTemperature(1.0);
        final ClusterChainGaussian refCluster = new ClusterChainGaussian(meanPosition, standardDeviation);
        refCluster.setTemperature(1.0);
        final double refIntegral = 1.0;
        System.out.println("reference integral: " + refIntegral);

        long blockSize = steps / 1000;
        System.out.println(steps + " steps (" + (steps / blockSize) + " blocks of " + blockSize + ")");
        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space, new SpeciesSpheresMono(space, new ElementSimple("A")), nPoints, 1, refCluster, targetCluster);
//        sim.setRandom(new RandomMersenneTwister(2));
//        int[] seeds = sim.getRandomSeeds();
//        System.out.println("Random Seeds:");
//        for (int i = 0; i < seeds.length; i++) {
//            System.out.println(seeds[i]);
//        }
        sim.setBoxLengths(params.targetL, params.targetL);
        sim.init();
        double factor = sim.boxLengths[1]/(sigmaB*nBig + sigmaS*nSmall);
//        System.out.println("Factor: " + factor);
        meanPosition[0][0] = -sim.boxLengths[1] / 2 + factor * 0.5 * pairSigma[0][0];
//        System.out.println("0 " + meanPosition[0][0]);
        sim.box[1].getLeafList().get(0).getPosition().setX(0, meanPosition[0][0]);
//        System.out.println("P0 " + sim.box[1].getLeafList().get(0).getPosition().getX(0));
        for (int i = 1; i < sim.box[1].getLeafList().size(); i++) {
            meanPosition[i][0] = sim.box[1].getLeafList().get(i - 1).getPosition().getX(0) + factor * pairSigma[i][i - 1];
            sim.box[1].getLeafList().get(i).getPosition().setX(0, meanPosition[i][0]);
//            System.out.println(i + " " + meanPosition[i][0]);
        }

//        for (int i = 0; i < nPoints; i++){
//            for (int j = 0; j < D; j++){
//                meanPosition[i][j] = sim.box[1].getLeafList().get(i).getPosition().getX(j);
//            }
//        }
        sim.integratorOS.setDoAdjustOnTime(true);

        sim.integrators[0].getMoveManager().removeMCMove(sim.mcMoveTranslate[0]);
        MCMoveClusterGaussian mcMoveClusterGaussian = new MCMoveClusterGaussian(sim.getRandom(), space, meanPosition, standardDeviation);
        sim.integrators[0].getMoveManager().addMCMove(mcMoveClusterGaussian);
        mcMoveClusterGaussian.setDoImposePBC(false);
        sim.integrators[1].getMoveManager().removeMCMove(sim.mcMoveTranslate[1]);
        MCMoveClusterFixCOM mcMoveClusterFixCOM = new MCMoveClusterFixCOM(sim.getRandom(), space);
        sim.integrators[1].getMoveManager().addMCMove(mcMoveClusterFixCOM);

//        ((MCMoveClusterAtomMulti) sim.mcMoveTranslate[0]).setStartAtom(0);
//        ((MCMoveClusterAtomMulti) sim.mcMoveTranslate[0]).setDoImposePBC(true);
//        ((MCMoveClusterAtomMulti) sim.mcMoveTranslate[1]).setStartAtom(0);
//        ((MCMoveClusterAtomMulti) sim.mcMoveTranslate[1]).setDoImposePBC(true);
//        MCMoveClusterAtomMulti targetBigMove = new MCMoveClusterAtomMulti(sim.getRandom(), space);
//        targetBigMove.setStepSize(params.targetL / 2);
//        targetBigMove.setStepSizeMin(params.targetL / 2);
//        sim.integrators[1].getMoveManager().addMCMove(targetBigMove);
//        ((MCMoveStepTracker) sim.mcMoveTranslate[1].getTracker()).setAcceptanceTarget(targetAcceptance);

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
//                    if (idx < nBig) return 1;
//                    return q;
                    return sigma[idx];
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

        System.out.println("equilibration finished");

        if (refFrac >= 0) {
            sim.integratorOS.setRefStepFraction(refFrac);
            sim.integratorOS.setAdjustStepFraction(false);
        }

        sim.integratorOS.setNumSubSteps((int) blockSize);
        sim.setAccumulatorBlockSize(blockSize);
        sim.ai.setMaxSteps(steps / blockSize);
        System.out.println("MC Move step size (target) " + mcMoveClusterFixCOM.getStepSize());

        MeterPercolation meter = new MeterPercolation(sim.box[1]);
        AccumulatorAverageFixed accumulator = new AccumulatorAverageFixed(blockSize);
        sim.integrators[1].getEventManager().addListener(new DataPumpListener(meter, accumulator));

        sim.getController().actionPerformed();
        long t2 = System.currentTimeMillis();
        //System.out.println("MC Move target big step acceptance " + targetBigMove.getTracker().acceptanceProbability() + " (" + targetBigMove.getStepSize() + ")");

        System.out.println("final reference step fraction " + sim.integratorOS.getIdealRefStepFraction());
        System.out.println("actual reference step fraction " + sim.integratorOS.getRefStepFraction());
        System.out.println("time reference step fraction " + sim.integratorOS.getRefTimeFraction());

        sim.printResults(refIntegral, null);

        DataGroup allYourBase = (DataGroup) accumulator.getData();
        IData averageData = allYourBase.getData(accumulator.AVERAGE.index);
        IData errorData = allYourBase.getData(accumulator.ERROR.index);

//        System.out.println();

        for (int i=0; i<sim.box[1].getLeafList().size(); i++){
//            System.out.print(i);
            for (int j=0; j<sim.getSpace().D(); j++){
                double avg = Math.sqrt(averageData.getValue(i*sim.getSpace().D() + j));
                double err = (errorData.getValue(i*sim.getSpace().D() + j))/(2*avg);
//                double err = errorData.getValue(i*sim.getSpace().D() + j);
                System.out.print(String.format("  %20.15e %9.4e", avg, err));
//                System.out.print(String.format("%20.15e ", avg));
            }
            System.out.println();
        }

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
        public double targetAcceptance = 0.2;
        public boolean writeRefPref = false;
        public String order = "BBSSSS";
//        public char[] order = {'B', 'B', 'S', 'S', 'S', 'S', 'S'};
    }
}
