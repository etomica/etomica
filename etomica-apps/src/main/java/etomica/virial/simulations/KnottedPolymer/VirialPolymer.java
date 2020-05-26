/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations.KnottedPolymer;

import etomica.action.IAction;
import etomica.box.Box;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.meter.MeterRadiusGyration;
import etomica.data.types.DataDouble;
import etomica.data.types.DataGroup;
import etomica.graphics.*;
import etomica.integrator.IntegratorListenerAction;
import etomica.math.SpecialFunctions;
import etomica.molecule.IMoleculeList;
import etomica.potential.IPotential;
import etomica.potential.P2HardSphere;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.units.dimensions.Null;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.virial.*;
import etomica.virial.simulations.SimulationVirial;

import javax.swing.*;
import java.awt.*;

/**
 * Calculation for virial coefficients of hard spheres
 * Edited for calculating protein 1HRB
 * <p>
 * CHANGES: change for star-polymers. May 7, 2018
 */
public class VirialPolymer {

    public static void main(String[] args) {

        VirialHSParam params = new VirialHSParam();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        } else {
            params.nPoints = 2;
            params.numSteps = 1000000L;
            params.ref = VirialHSParam.TREE;
            params.chainFrac = 0.2;
            params.ringFrac = 0.7;
            params.nPtsTabulated = 0;
        }

        final int nPoints = params.nPoints;
        long steps = params.numSteps;
        final int ref = params.ref;

        Space space = Space3D.getInstance();
        SpeciesStarPolymer species = new SpeciesStarPolymer(space, "./resource/f5L40_00.xyz");

//        final double sigmaTranslate = 10.0;
        final double sigmaTranslate = species.maxDiameter;
        final int nPtsTabulated = params.nPtsTabulated;

        double refDiameter = 5.0;
        double vhs = (4.0 / 3.0) * Math.PI * refDiameter * refDiameter * refDiameter;
        System.out.println("HS singly-connected sampling B" + nPoints);
        System.out.println("Max Diameter for this polymer is " + sigmaTranslate);


//        SpeciesStarPolymer species = new SpeciesStarPolymer(space, "./resource/two.xyz");

//        final P2HardSphere potentialHS = new P2HardSphere(space);
//        final PotentialGroup potentialGroup = new PotentialGroup(2, space);
//        potentialGroup.addPotential(potentialHS, new AtomType[]{species.getAtomType(0), species.getAtomType(0)});
//        MayerGeneral fRef = new MayerGeneral(potentialGroup);
        final P2HSPolymer potentialHS = new P2HSPolymer(2, space, new P2HardSphere(space));
        MayerGeneral fRef = new MayerGeneral(potentialHS);
        MayerFunction fRefPos = new MayerFunction() {
            public void setBox(Box box) {
            }

            public IPotential getPotential() {
                return null;
            }

            public double f(IMoleculeList pair, double r2, double beta) {
                return r2 < sigmaTranslate * sigmaTranslate ? 1 : 0;
            }
        };

        ClusterAbstract targetCluster = nPtsTabulated > 0 ? new ClusterWheatleyPartitionScreening(nPoints, fRef, nPtsTabulated) : new ClusterWheatleyHS(nPoints, fRef);

        targetCluster.setTemperature(1.0);

        ClusterAbstract refCluster = null;
        long numDiagrams;

        double ri = 0;
        if (ref == VirialHSParam.TREE) {
            System.out.println("using a tree reference");
            refCluster = new ClusterSinglyConnected(nPoints, fRefPos);
            numDiagrams = ((ClusterSinglyConnected) refCluster).numDiagrams();
            ri = numDiagrams * Math.pow(vhs, nPoints - 1);
        } else if (ref == VirialHSParam.CHAINS) {
            System.out.println("using a chain reference");
            refCluster = new ClusterChainHS(nPoints, fRefPos);
            numDiagrams = ((ClusterChainHS) refCluster).numDiagrams();
            ri = numDiagrams * Math.pow(vhs, nPoints - 1);
        }

        // (nPoints-1)! is simply not included by ClusterWheatley, so do that here.
        final double refIntegral = ri / SpecialFunctions.factorial(nPoints);
        System.out.println("reference integral: " + refIntegral);
        refCluster.setTemperature(1.0);

        ClusterAbstract[] targetDiagrams = new ClusterAbstract[]{targetCluster};

        System.out.println(steps + " steps");

        final SimulationVirial sim = new SimulationVirial(space, species,
                1.0, ClusterWeightAbs.makeWeightCluster(refCluster), refCluster, targetDiagrams);

//        final SimulationVirialOverlap sim = new SimulationVirialOverlap(space, new SpeciesFactoryStarPolymer(), 0, refCluster, targetCluster);


        sim.box.getBoundary().setBoxSize(Vector.of(new double[]{50, 50, 50}));
//        sim.box.getBoundary().setBoxSize(Vector.of(new double[]{2*sigmaTranslate, 2*sigmaTranslate, 2*sigmaTranslate}));
        MeterVirialBD meter = new MeterVirialBD(sim.allValueClusters);
        meter.setBox(sim.box);
        sim.setMeter(meter);
        AccumulatorAverageFixed accumulator = new AccumulatorAverageFixed(1000);
        sim.setAccumulator(accumulator);
        accumulator.setPushInterval(100000000);

        if (false) {

            sim.box.getBoundary().setBoxSize(Vector.of(new double[]{100, 100, 100}));
//            sim.box.getBoundary().setBoxSize(Vector.of(new double[]{2*sigmaTranslate, 2*sigmaTranslate, 2*sigmaTranslate}));
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
            DisplayBox displayBox0 = simGraphic.getDisplayBox(sim.box);
//            displayBox0.setPixelUnit(new Pixel(300.0/size));
            displayBox0.setShowBoundary(false);
            ((DisplayBoxCanvasG3DSys) displayBox0.canvas).setBackgroundColor(Color.WHITE);

            ColorScheme colorScheme = new ColorSchemeRandomByMolecule(sim, sim.box, sim.getRandom());
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
            IDataInfo dataInfoSteps = new DataDouble.DataInfoDouble("B" + nPoints, Null.DIMENSION);
            stepsBox.putDataInfo(dataInfoSteps);
            IDataInfo dataInfo = new DataDouble.DataInfoDouble("B" + nPoints, Null.DIMENSION);
            averageBox.putDataInfo(dataInfo);
            averageBox.setLabel("average");
            errorBox.putDataInfo(dataInfo);
            errorBox.setLabel("error");
            errorBox.setPrecision(2);
            sim.integrator.getEventManager().addListener(new IntegratorListenerAction(pushAnswer, 1000));
            simGraphic.makeAndDisplayFrame();

            return;
        }

        long t1 = System.currentTimeMillis();

        sim.ai.setMaxSteps(steps);
        sim.getController().actionPerformed();
        long t2 = System.currentTimeMillis();

//        if (!Double.isNaN(litHSB)) System.out.println("lit value "+litHSB);

        DataGroup allYourBase = (DataGroup) accumulator.getData();
        IData averageData = allYourBase.getData(accumulator.AVERAGE.index);
        IData errorData = allYourBase.getData(accumulator.ERROR.index);

        double AcceptProbTrans = sim.mcMoveTranslate.getTracker().acceptanceProbability();
        double AcceptProbRotate = sim.mcMoveRotate.getTracker().acceptanceProbability();
        System.out.println();
        System.out.print(String.format("Acceptance Probability of McMoveTranslate: %f\n", AcceptProbTrans));
        System.out.print(String.format("Acceptance Probability of McMoveRotate: %f\n", AcceptProbRotate));
        System.out.println();

        double avg = averageData.getValue(0);
        double err = errorData.getValue(0);
        System.out.print(String.format("target average: %20.15e error: %9.4e\n", avg, err));

        System.out.println();

        System.out.print(String.format("abs average: %20.15e  error: %9.4e\n", avg * refIntegral, err * Math.abs(refIntegral)));

//        sim.printResults(refIntegral);
        MeterRadiusGyration meterGyration = new MeterRadiusGyration(sim.getSpace());
        meterGyration.setBox(sim.box);
        sim.setMeter(meterGyration);
        System.out.println(String.format("Radius of Gyration: %.8f\n", sim.meter.getData().getValue(0)));

        System.out.println("time: " + (t2 - t1) / 1000.0);


        if (targetCluster instanceof ClusterWheatleyPartitionScreening) {
            ClusterWheatleyPartitionScreening tc = (ClusterWheatleyPartitionScreening) targetCluster;
            if (tc.sigCounter != null) {
                for (int i = 3; i < 6; i++) {
                    System.out.println("size " + i);
                    tc.sigCounter[i - 1].print();
                }
                tc.sigCounter[nPoints - 1].print();
                System.out.println("Fraction tabulated for each n");
                for (int i = 0; i < nPoints; i++) {
                    System.out.println(tc.sigCounter[i].getn() + "\t" + tc.sigCounter[i].fractionTabulated() + "\t" + tc.sigCounter[i].getEntries());
                }
            }
        }
    }

    /**
     * Inner class for parameters
     */
    public static class VirialHSParam extends ParameterBase {
        public int nPoints = 10;
        public long numSteps = 100000000;
        public static final int TREE = 0, CHAINS = 1;
        public int ref = TREE;
        public double chainFrac = 1;
        public double ringFrac = 0.5;
        public int nPtsTabulated = 0;
    }

}
