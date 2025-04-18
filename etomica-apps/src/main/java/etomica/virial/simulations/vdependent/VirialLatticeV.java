/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations.vdependent;

import etomica.action.IAction;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.chem.elements.ElementSimple;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.types.DataDouble;
import etomica.data.types.DataGroup;
import etomica.graphics.*;
import etomica.integrator.IntegratorListenerAction;
import etomica.molecule.IMoleculeList;
import etomica.potential.P2SquareWell;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularNonperiodic;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.species.ISpecies;
import etomica.species.SpeciesGeneral;
import etomica.units.dimensions.Null;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.virial.MayerFunction;
import etomica.virial.MayerGeneralSpherical;
import etomica.virial.MeterVirial;
import etomica.virial.cluster.*;
import etomica.virial.mcmove.MCMoveClusterAtomHSChain;
import etomica.virial.mcmove.MCMoveClusterAtomHSTree;
import etomica.virial.mcmove.MCMoveClusterAtomInBox;
import etomica.virial.mcmove.RandomPositionLattice;
import etomica.virial.simulations.SimulationVirial;
import etomica.virial.wheatley.ClusterWheatleySoft;

import javax.swing.*;
import java.awt.*;

/**
 * Calculation for virial coefficients of hard sphere mixture
 *
 * @author Arpit, Pavan and Andrew
 */
public class VirialLatticeV {

    public static void main(String[] args) {

        VirialHSParam params = new VirialHSParam();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        } else {
            params.nPoints = 2;
            params.numSteps = 10000000L;
            params.ref = ReferenceChoice.CHAIN_TREE;
            params.chainFrac = 0.5;
            params.D = 2;
            params.L = 2;
        }
        final int nPoints = params.nPoints;
        long steps = params.numSteps;
        final ReferenceChoice ref = params.ref;
        final double chainFrac = params.chainFrac;
        final int D = params.D;
        final double L = params.L;
        final double temperature = params.T;
        System.out.println("Number of points: " + nPoints);
        System.out.println("Dimensions: " + D);
        System.out.println("Box Length: " + L);
        Space space = Space.getInstance(D);
        long t1 = System.currentTimeMillis();

        System.out.println("B" + nPoints);

        Boundary b = (L < Double.POSITIVE_INFINITY && L > 0) ? new BoundaryRectangularPeriodic(space, L) : new BoundaryRectangularNonperiodic(space);

        double epsilon = L > 2 ? 1.0 : 2.0;
        MayerVDependent fTarget = new MayerVDependent(new MayerGeneralSpherical(P2SquareWell.makePotential(0.9, 1.2/0.9, epsilon)), b);
        MayerVDependent fRefPos = new MayerVDependent(new MayerFunction() {
            @Override
            public double f(IMoleculeList pair, double r2, double beta) {
                return r2 < 1.2 ? 1 : 0;
            }

            @Override
            public void setBox(Box box) {}
        }, b);

        final ClusterAbstract targetCluster = new ClusterWheatleySoft(nPoints, fTarget, 1e-12);
        targetCluster.setTemperature(temperature);
        final ClusterWeightUmbrella refCluster;
        if (ref == ReferenceChoice.TREE) {
            System.out.println("using a tree reference");
            ClusterSinglyConnected ct = new ClusterSinglyConnected(nPoints, fRefPos);
            refCluster = new ClusterWeightUmbrella(new ClusterAbstract[]{ct});
            refCluster.setWeightCoefficients(new double[]{1.0 / ct.numDiagrams()});
        } else if (ref == ReferenceChoice.CHAINS) {
            System.out.println("using a chain reference");
            ClusterChainHS cc = new ClusterChainHS(nPoints, fRefPos);
            refCluster = new ClusterWeightUmbrella(new ClusterAbstract[]{cc});
            refCluster.setWeightCoefficients(new double[]{1.0 / cc.numDiagrams()});
        } else if (ref == ReferenceChoice.CHAIN_TREE) {
            System.out.println("using a chain/tree reference (" + chainFrac + " chains)");
            ClusterChainHS cc = new ClusterChainHS(nPoints, fRefPos);
            ClusterSinglyConnected ct = new ClusterSinglyConnected(nPoints, fRefPos);
            refCluster = new ClusterWeightUmbrella(new ClusterAbstract[]{cc, ct});
            refCluster.setWeightCoefficients(new double[]{chainFrac / cc.numDiagrams(), (1 - chainFrac) / ct.numDiagrams()});
        } else if (ref == ReferenceChoice.RANDOM) {
            System.out.println("using random particles in box reference");
            ClusterConstant cc = new ClusterConstant(nPoints, Math.pow(1.0/L, 3.0*(nPoints-1)));
            refCluster = new ClusterWeightUmbrella(new ClusterAbstract[]{cc});
            refCluster.setWeightCoefficients(new double[]{1.0});
        } else {
            throw new RuntimeException();
        }

        refCluster.setTemperature(1.0);

        ClusterAbstract[] targetDiagrams = new ClusterAbstract[]{targetCluster};

        System.out.println(steps + " steps ");
        final SimulationVirial sim = new SimulationVirial(space, new ISpecies[]{SpeciesGeneral.monatomic(space, AtomType.element(new ElementSimple("A")))}, new int[]{nPoints}, 1.0, ClusterWeightAbs.makeWeightCluster(refCluster), refCluster, targetDiagrams);
        if (L > 0 && L < Double.POSITIVE_INFINITY) sim.setBoxLength(L);
        sim.init();
        MeterVirial meter = new MeterVirial(sim.allValueClusters);
        meter.setBox(sim.box);
        sim.setMeter(meter);
        AccumulatorAverageFixed accumulator = new AccumulatorAverageFixed(1000);
        sim.setAccumulator(accumulator);
        accumulator.setPushInterval(1000000);

        sim.integrator.getMoveManager().removeMCMove(sim.mcMoveTranslate);
        RandomPositionLattice positionSource = new RandomPositionLattice(space, sim.getRandom(), L, 1);
        if (ref == ReferenceChoice.TREE) {
            MCMoveClusterAtomHSTree mcMoveHS = new MCMoveClusterAtomHSTree(sim.getRandom(), sim.box, 1);
            mcMoveHS.setPositionSource(positionSource);
            sim.integrator.getMoveManager().addMCMove(mcMoveHS);
        } else if (ref == ReferenceChoice.CHAINS) {
            MCMoveClusterAtomHSChain mcMoveHS = new MCMoveClusterAtomHSChain(sim.getRandom(), sim.box, 1);
            mcMoveHS.setPositionSource(positionSource);
            sim.integrator.getMoveManager().addMCMove(mcMoveHS);
        } else if (ref == ReferenceChoice.CHAIN_TREE) {
            MCMoveClusterAtomHSTree mcMoveHST = new MCMoveClusterAtomHSTree(sim.getRandom(), sim.box, 1);
            mcMoveHST.setPositionSource(positionSource);
            sim.integrator.getMoveManager().addMCMove(mcMoveHST);
            sim.integrator.getMoveManager().setFrequency(mcMoveHST, 1 - chainFrac);
            MCMoveClusterAtomHSChain mcMoveHSC = new MCMoveClusterAtomHSChain(sim.getRandom(), sim.box, 1);
            mcMoveHSC.setPositionSource(positionSource);
            sim.integrator.getMoveManager().addMCMove(mcMoveHSC);
            sim.integrator.getMoveManager().setFrequency(mcMoveHSC, chainFrac);
        } else if (ref == ReferenceChoice.RANDOM) {
            MCMoveClusterAtomInBox mcMoveHS = new MCMoveClusterAtomInBox(sim.getRandom(), sim.box);
            sim.integrator.getMoveManager().addMCMove(mcMoveHS);
        }
        else {
            throw new RuntimeException("Unknown reference");
        }

        if (false) {
            if (L == Double.POSITIVE_INFINITY || L == 0)
                sim.box.getBoundary().setBoxSize(space.makeVector(new double[]{10, 10, 10}));
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
            DisplayBox displayBox0 = simGraphic.getDisplayBox(sim.box);
//            displayBox0.setPixelUnit(new Pixel(300.0/size));
            if (L == Double.POSITIVE_INFINITY || L == 0) {
                displayBox0.setShowBoundary(false);
                ((DisplayBoxCanvasG3DSys) displayBox0.canvas).setBackgroundColor(Color.WHITE);
            }

            simGraphic.makeAndDisplayFrame();

            accumulator.setBlockSize(1000);

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

            IAction pushAnswer = new IAction() {
                public void actionPerformed() {
                    IData avgData = accumulator.getData();
                    double avg = ((DataGroup) avgData).getData(accumulator.AVERAGE.index).getValue(0);
                    double error = ((DataGroup) avgData).getData(accumulator.ERROR.index).getValue(0);
                    data.x = avg;
                    averageBox.putData(data);
                    data.x = error;
                    errorBox.putData(data);

                    data.x = sim.integrator.getStepCount();
                    stepsBox.putData(data);
                }

                DataDouble data = new DataDouble();
            };
            IDataInfo dataInfoSteps = new DataDouble.DataInfoDouble("B" + nPoints, Null.DIMENSION);
            stepsBox.putDataInfo(dataInfoSteps);
            IDataInfo dataInfo = new DataDouble.DataInfoDouble("B" + nPoints, Null.DIMENSION);
            averageBox.putDataInfo(dataInfo);
            averageBox.setLabel("average");
            errorBox.putDataInfo(dataInfo);
            errorBox.setLabel("error");
            errorBox.setPrecision(2);
            sim.integrator.getEventManager().addListener(new IntegratorListenerAction(pushAnswer, 1000));

            return;
        }


        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, steps));


        DataGroup allYourBase = (DataGroup) accumulator.getData();
        IData averageData = allYourBase.getData(accumulator.AVERAGE.index);
        IData errorData = allYourBase.getData(accumulator.ERROR.index);

        System.out.println();

        double avg = averageData.getValue(1);
        double err = errorData.getValue(1);

        System.out.printf("target average: %20.15e error: %9.4e\n", avg, err);

        long t2 = System.currentTimeMillis();
        System.out.println("time: " + (t2 - t1) / 1000.0);
    }

    public enum ReferenceChoice {
        TREE, CHAINS, CHAIN_TREE, RANDOM
    }

    /**
     * Inner class for parameters
     */
    public static class VirialHSParam extends ParameterBase {
        public int nPoints = 3;
        public long numSteps = 100000000;
        public ReferenceChoice ref = ReferenceChoice.CHAIN_TREE;
        public double chainFrac = 0.5;
        public int D = 3;
        // set L to 0 or infinity to get a coefficient without PBC
        public int L = 0;
        public double T = 1;
    }

}
