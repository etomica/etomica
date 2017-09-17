package etomica.virial.simulations;

import etomica.action.IAction;
import etomica.atom.IAtom;
import etomica.chem.elements.ElementSimple;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.IData;
import etomica.data.IEtomicaDataInfo;
import etomica.data.types.DataDouble;
import etomica.data.types.DataGroup;
import etomica.graphics.*;
import etomica.listener.IntegratorListenerAction;
import etomica.space.Space;
import etomica.species.SpeciesSpheresMono;
import etomica.units.dimensions.Null;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.virial.*;

import javax.swing.*;
import java.awt.*;

/**
 * Calculation for virial coefficients of hard sphere mixture
 *
 * @author Pavan and Andrew
 */
public class VirialHSMixture {

    public static void main(String[] args) {

        VirialHSParam params = new VirialHSParam();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        } else {
            params.nPoints = 5;
            params.numSteps = 10000000L;
            params.ref = VirialHSParam.CHAIN_TREE;
            params.chainFrac = 0.5;
            params.q = 0.5;
            params.nonAdd = 0.05;
            params.cff = 1;
            params.D = 2;
        }
        final int nSmall = params.cff;
        final int nPoints = params.nPoints;
        long steps = params.numSteps;
        final double q = params.q;
        final int ref = params.ref;
        System.out.println("   q: " + q);
        final double sigmaB = 1;
        final double sigmaS = q * sigmaB;
        final double[][] nonAdd = new double[2][2];
        nonAdd[0][1] = nonAdd[1][0] = params.nonAdd;
        final double chainFrac = params.chainFrac;
        final int D = params.D;
        System.out.println("Number of points: " + nPoints);
        System.out.println("Dimensions: " + D);
        System.out.println("nonAdd: " + params.nonAdd);
        Space space = Space.getInstance(D);
        long t1 = System.currentTimeMillis();

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
                        double x = nonAdd == null ? 0 : nonAdd[i][j];
                        pairSigma[jjj][iii] = pairSigma[iii][jjj] = (1 + x) * (sigma[i] + sigma[j]) / 2;
                        jjj++;
                    }
                }
                iii++;
            }
        }

        MayerHSMixture fTarget = MayerHSMixture.makeTargetF(nPoints, pairSigma);
        MayerFunction fRefPos = MayerHSMixture.makeReferenceF(space, nPoints, pairSigma);

        final ClusterAbstract targetCluster = new ClusterWheatleyHS(nPoints, fTarget);
        targetCluster.setTemperature(1.0);
        final ClusterAbstract refCluster;
        double ri;
        if (ref == VirialHSParam.TREE) {
            System.out.println("using a tree reference");
            ClusterSinglyConnected ct = new ClusterSinglyConnected(nPoints, fRefPos);
            refCluster = new ClusterWeightUmbrella(new ClusterAbstract[]{ct});
            ((ClusterWeightUmbrella) refCluster).setWeightCoefficients(new double[]{1.0 / ct.numDiagrams()});
            ri = 1;
        } else if (ref == VirialHSParam.CHAINS) {
            System.out.println("using a chain reference");
            ClusterChainHS cc = new ClusterChainHS(nPoints, fRefPos);
            refCluster = new ClusterWeightUmbrella(new ClusterAbstract[]{cc});
            ((ClusterWeightUmbrella) refCluster).setWeightCoefficients(new double[]{1.0 / cc.numDiagrams()});
            ri = 1;
        } else if (ref == VirialHSParam.CHAIN_TREE) {
            System.out.println("using a chain/tree reference (" + chainFrac + " chains)");
            ClusterChainHS cc = new ClusterChainHS(nPoints, fRefPos);
            ClusterSinglyConnected ct = new ClusterSinglyConnected(nPoints, fRefPos);
            refCluster = new ClusterWeightUmbrella(new ClusterAbstract[]{cc, ct});
            ((ClusterWeightUmbrella) refCluster).setWeightCoefficients(new double[]{chainFrac / cc.numDiagrams(), (1 - chainFrac) / ct.numDiagrams()});
            ri = 1;
        } else {
            throw new RuntimeException();
        }

        final double refIntegral = ri;
        System.out.println("reference integral: " + refIntegral);
        refCluster.setTemperature(1.0);

        ClusterAbstract[] targetDiagrams = new ClusterAbstract[]{targetCluster};

        System.out.println(steps + " steps ");
        final SimulationVirial sim = new SimulationVirial(space, new SpeciesSpheresMono(space, new ElementSimple("A")), 1.0, ClusterWeightAbs.makeWeightCluster(refCluster), refCluster, targetDiagrams, false);
        MeterVirialBD meter = new MeterVirialBD(sim.allValueClusters);
        meter.setBox(sim.box);
        sim.setMeter(meter);
        AccumulatorAverageFixed accumulator = new AccumulatorAverageFixed(1000);
        sim.setAccumulator(accumulator);
        accumulator.setPushInterval(1000000);

        sim.integrator.getMoveManager().removeMCMove(sim.mcMoveTranslate);
        if (ref == VirialHSParam.TREE) {
            MCMoveClusterAtomHSTreeMix mcMoveHS = new MCMoveClusterAtomHSTreeMix(sim.getRandom(), space, pairSigma);
            sim.integrator.getMoveManager().addMCMove(mcMoveHS);
        } else if (ref == VirialHSParam.CHAINS) {
            MCMoveClusterAtomHSChainMix mcMoveHS = new MCMoveClusterAtomHSChainMix(sim.getRandom(), space, pairSigma);
            sim.integrator.getMoveManager().addMCMove(mcMoveHS);
        } else if (ref == VirialHSParam.CHAIN_TREE) {
            MCMoveClusterAtomHSTreeMix mcMoveHST = new MCMoveClusterAtomHSTreeMix(sim.getRandom(), space, pairSigma);
            sim.integrator.getMoveManager().addMCMove(mcMoveHST);
            sim.integrator.getMoveManager().setFrequency(mcMoveHST, 1 - chainFrac);
            MCMoveClusterAtomHSChainMix mcMoveHSC = new MCMoveClusterAtomHSChainMix(sim.getRandom(), space, pairSigma);
            sim.integrator.getMoveManager().addMCMove(mcMoveHSC);
            sim.integrator.getMoveManager().setFrequency(mcMoveHSC, chainFrac);
        }

        if (false) {
            sim.box.getBoundary().setBoxSize(space.makeVector(new double[]{10, 10, 10}));
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, space, sim.getController());
            DisplayBox displayBox0 = simGraphic.getDisplayBox(sim.box);
//            displayBox0.setPixelUnit(new Pixel(300.0/size));
            displayBox0.setShowBoundary(false);
            ((DisplayBoxCanvasG3DSys) displayBox0.canvas).setBackgroundColor(Color.WHITE);


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

            final JPanel panelParentGroup = new JPanel(new java.awt.GridBagLayout());
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

            IAction pushAnswer = new IAction() {
                public void actionPerformed() {
                    IData avgData = sim.accumulator.getData();
                    double ratio = ((DataGroup) avgData).getData(sim.accumulator.RATIO.index).getValue(1);
                    double error = ((DataGroup) avgData).getData(sim.accumulator.RATIO_ERROR.index).getValue(1);
                    data.x = ratio * refIntegral;
                    averageBox.putData(data);
                    data.x = error * Math.abs(refIntegral);
                    errorBox.putData(data);

                    data.x = sim.integrator.getStepCount();
                    stepsBox.putData(data);
                }

                DataDouble data = new DataDouble();
            };
            IEtomicaDataInfo dataInfoSteps = new DataDouble.DataInfoDouble("B" + nPoints, Null.DIMENSION);
            stepsBox.putDataInfo(dataInfoSteps);
            IEtomicaDataInfo dataInfo = new DataDouble.DataInfoDouble("B" + nPoints, Null.DIMENSION);
            averageBox.putDataInfo(dataInfo);
            averageBox.setLabel("average");
            errorBox.putDataInfo(dataInfo);
            errorBox.setLabel("error");
            errorBox.setPrecision(2);
            sim.integrator.getEventManager().addListener(new IntegratorListenerAction(pushAnswer, 1000));

            return;
        }


        sim.ai.setMaxSteps(steps);
        sim.getController().actionPerformed();


        DataGroup allYourBase = (DataGroup) accumulator.getData();
        IData averageData = allYourBase.getData(accumulator.AVERAGE.index);
        IData errorData = allYourBase.getData(accumulator.ERROR.index);

        System.out.println();

        double avg = averageData.getValue(0);
        double err = errorData.getValue(0);


        System.out.print(String.format("target average: %20.15e error: %9.4e\n", avg, err));

        System.out.println();

        System.out.println(String.format("abs average: %20.15e  error: %9.4e\n", avg * refIntegral, err * Math.abs(refIntegral)));

        long t2 = System.currentTimeMillis();
        System.out.println("time:" + (t2 - t1) / 1000.0);
        System.out.println("#################################");
    }


    /**
     * Inner class for parameters
     */
    public static class VirialHSParam extends ParameterBase {
        public int nPoints = 3;
        public long numSteps = 100000000;
        public static final int TREE = 0, CHAINS = 1, CHAIN_TAIL = 4, CHAIN_TREE = 5, CRINGS = 6, RING_TREE = 7, RINGS = 8, RING_CHAIN_TREES = 9;
        public int ref = CHAIN_TREE;
        public double chainFrac = 0.5;
        public double q = 0.1;
        public int cff = 0;
        public double nonAdd = 0;
        public int D = 3;
    }

}
