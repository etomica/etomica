/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations.theta;

import etomica.action.IAction;
import etomica.action.activity.ActivityIntegrate;
import etomica.action.controller.Activity;
import etomica.action.controller.Controller;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.data.DataPumpListener;
import etomica.data.DataSourceCountSteps;
import etomica.data.IDataInfo;
import etomica.data.types.DataDouble;
import etomica.graphics.*;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.mcmove.MCMoveStepTracker;
import etomica.math.SpecialFunctions;
import etomica.molecule.IMoleculeList;
import etomica.potential.IPotential2;
import etomica.potential.P2LennardJones;
import etomica.potential.PotentialMasterBonding;
import etomica.potential.PotentialMoleculePair;
import etomica.potential.compute.PotentialCompute;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesBuilder;
import etomica.species.SpeciesManager;
import etomica.starpolymer.ConformationStarPolymerGraft;
import etomica.units.Pixel;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.*;
import etomica.util.Constants.CompassDirection;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.util.collections.IntArrayList;
import etomica.virial.MayerFunction;
import etomica.virial.MayerGeneral;
import etomica.virial.cluster.*;
import etomica.virial.mcmove.*;
import etomica.virial.simulations.SimulationVirialOverlap2;

import javax.swing.*;
import java.awt.*;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Mayer sampling simulation for alkanes using the TraPPE-United Atoms force field.
 *   M.G. Martin and J.I. Siepmann, "Transferable Potentials for Phase
 *   Equilibria. 1. United-Atom Description of n-Alkanes," J. Phys. Chem. B
 *   102, 2569-2577 (1998)
 *
 *   modified from VirialAlkaneFlex2 so that eovererr can be used
 *   March 2013
 */
public class VirialStar {

    public static ClusterSum makeB2Cluster(MayerFunction f) {
        return new ClusterSum(new ClusterBonds[]{new ClusterBonds(2, new int[][][]{{{0,1}}})}, new double[]{1}, new MayerFunction[]{f});
    }

    public static void main(String[] args) {
        VirialStarParams params = new VirialStarParams();
        boolean isCommandline = args.length > 0;
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        } else {
            params.nPoints = 2;
            params.armLength = 1;
            params.numArms = 4;
            params.temperature = 4.6;
            params.numSteps = 10000000;
        }
        final int nPoints = params.nPoints;
        int numArms = params.numArms;
        int armLength = params.armLength;
        double temperature = params.temperature;
        long steps = params.numSteps;
        double refFreq = params.refFreq;
        double sigmaHSRef = 1.5 + 0.15*(2*armLength);

        Space space = Space3D.getInstance();

        ConformationStarPolymerGraft conf = new ConformationStarPolymerGraft(Space3D.getInstance(), numArms, armLength);
        conf.setSigma0(1);
        AtomType type = AtomType.simple("A");
        ISpecies species = new SpeciesBuilder(Space3D.getInstance())
                .addCount(type, 1+numArms*armLength)
                .withConformation(conf).build();
        SpeciesManager sm = new SpeciesManager.Builder().addSpecies(species).build();

        PotentialMoleculePair pTarget = new PotentialMoleculePair(space, sm);
        System.out.println(numArms+" arms of length "+armLength+" B"+nPoints+" at T = "+temperature);
        IPotential2 p2 = new P2LennardJones(1, 1);

        MayerGeneral fTarget = new MayerGeneral(pTarget);
        MayerFunction fRefPos = new MayerFunction() {
            public void setBox(Box box) {
            }

            public double f(IMoleculeList pair, double r2, double beta) {
                return r2 < sigmaHSRef * sigmaHSRef ? 1 : 0;
            }
        };

        boolean flex = nPoints > 2;
        VirialDiagrams alkaneDiagrams = new VirialDiagrams(nPoints, false, flex);
        alkaneDiagrams.setDoReeHoover(false);
        ClusterAbstract targetCluster = null;

        targetCluster = alkaneDiagrams.makeVirialCluster(fTarget);
        ClusterChainHS refCluster = new ClusterChainHS(nPoints, fRefPos);
        double vhs = (4.0 / 3.0) * Math.PI * sigmaHSRef * sigmaHSRef * sigmaHSRef;
        final double refIntegral = SpecialFunctions.factorial(nPoints) / 2 * Math.pow(vhs, nPoints - 1);

        targetCluster.setTemperature(temperature);
        refCluster.setTemperature(temperature);

        System.out.println("sigmaHSRef: "+sigmaHSRef);
        // eovererr expects this string, BnHS
        System.out.println("B"+nPoints+"HS: "+refIntegral);

        PotentialMasterBonding.FullBondingInfo bondingInfo = new PotentialMasterBonding.FullBondingInfo(sm) {
            @Override
            public boolean skipBondedPair(boolean isPureAtoms, IAtom iAtom, IAtom jAtom) {
                return false;
            }
        };

        // we need to do this to convince the system that the molecules are not rigid
        // if bondingInfo thinks molecules are rigid then intramolecular LJ will not be computed
        IPotential2 pBonding = new IPotential2() {
            @Override
            public double getRange() { return 2; }
            @Override
            public void u012add(double r2, double[] u012) { }
        };
        List<int[]> pairs = new ArrayList<>();
        for (int i=0; i<numArms; i++) {
            pairs.add(new int[]{0,1+armLength*i});
            for (int j = 0; j < armLength-1; j++) {
                pairs.add(new int[]{1+armLength*i+j, 1+armLength*i+j+1});
            }
        }
        bondingInfo.setBondingPotentialPair(species, pBonding, pairs);

        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space, sm,
                new int[]{flex ? (nPoints+1) : nPoints},temperature, refCluster, targetCluster);
        sim.setDoWiggle(false);
        sim.setBondingInfo(bondingInfo);
        sim.setIntraPairPotentials(pTarget.getAtomPotentials());
        sim.init();
        System.out.println("random seeds: "+ Arrays.toString(sim.getRandomSeeds()));

        PotentialMasterBonding.FullBondingInfo bondingInfodk = new PotentialMasterBonding.FullBondingInfo(sm);

        sim.integrators[0].getMoveManager().removeMCMove(sim.mcMoveTranslate[0]);
        MCMoveClusterMoleculeHSChain mcMoveHSC = new MCMoveClusterMoleculeHSChain(sim.getRandom(), sim.box[0], sigmaHSRef);
        sim.integrators[0].getMoveManager().addMCMove(mcMoveHSC);
        sim.accumulators[0].setBlockSize(1);

        int[] constraintMap = new int[nPoints+1];
        for (int i=0; i<nPoints; i++) {
            constraintMap[i] = i;
        }
        if (flex) {
            constraintMap[nPoints] = 0;
            mcMoveHSC.setConstraintMap(constraintMap);
            ((MCMoveClusterMoleculeMulti)sim.mcMoveTranslate[1]).setConstraintMap(constraintMap);
            ((MCMoveClusterRotateMoleculeMulti)sim.mcMoveRotate[0]).setConstraintMap(constraintMap);
            ((MCMoveClusterRotateMoleculeMulti)sim.mcMoveRotate[1]).setConstraintMap(constraintMap);
        }


        if (refFreq >= 0) {
            sim.integratorOS.setAdjustStepFraction(false);
            sim.integratorOS.setRefStepFraction(refFreq);
        }
      
        sim.integratorOS.setNumSubSteps(1000);
        sim.integratorOS.setAggressiveAdjustStepFraction(true);
        System.out.println(steps+" steps (1000 blocks of "+steps/1000+")");
        steps /= 1000;

        pTarget.setAtomPotential(type, type, p2);

        sim.integratorOS.setNumSubSteps(1000);

        MCMoveClusterAngle[] angleMoves = null;
        MCMoveClusterStretch[] stretchMoves = null;

        IntArrayList[] bonding = new IntArrayList[1+numArms*armLength];
        bonding[0] = new IntArrayList(numArms);

        for (int i=0; i<numArms; i++) {
            bonding[0].add(1+i*armLength);
            if (armLength > 1) {
                bonding[1 + i * armLength] = new IntArrayList(new int[]{0, 1 + i * armLength + 1});
            }
            else {
                bonding[1 + i * armLength] = new IntArrayList(new int[]{0});
            }
            for (int j = 1; j < armLength-1; j++) {
                bonding[1+armLength*i+j] = new IntArrayList(new int []{1+armLength*i+j-1, 1+armLength*i+j+1});
            }
            if (armLength > 1) {
                bonding[1 + (i+1) * armLength - 1] = new IntArrayList(new int[]{1 + (i+1) * armLength - 2});
            }
        }

        PotentialCompute pc0 = sim.integrators[0].getPotentialCompute();
        PotentialCompute pc1 = sim.integrators[1].getPotentialCompute();

        angleMoves = new MCMoveClusterAngle[2];
        angleMoves[0] = new MCMoveClusterAngle(pc0, space, bonding, sim.getRandom(), 1);
        angleMoves[0].setBox(sim.box[0]);
        angleMoves[0].setConstraintMap(constraintMap);
        sim.integrators[0].getMoveManager().addMCMove(angleMoves[0]);
        angleMoves[1] = new MCMoveClusterAngle(pc1, space, bonding, sim.getRandom(), 1);
        angleMoves[1].setBox(sim.box[1]);
        angleMoves[1].setConstraintMap(constraintMap);
        sim.integrators[1].getMoveManager().addMCMove(angleMoves[1]);

        MCMoveClusterShuffle shuffleMove0 = null, shuffleMove1 = null;
        if (armLength > 5) {
            shuffleMove0 = new MCMoveClusterShuffle(pc0, space, sim.getRandom());
            shuffleMove0.setBox(sim.box[0]);
            shuffleMove0.setBonding(bonding);
            shuffleMove0.setStepSizeMax(armLength - 2);
            shuffleMove0.setConstraintMap(constraintMap);
            sim.integrators[0].getMoveManager().addMCMove(shuffleMove0);
            ((MCMoveStepTracker) shuffleMove0.getTracker()).setAcceptanceTarget(0.3);
            shuffleMove1 = new MCMoveClusterShuffle(pc1, space, sim.getRandom());
            shuffleMove1.setBonding(bonding);
            shuffleMove1.setBox(sim.box[1]);
            shuffleMove1.setStepSizeMax(armLength - 2);
            shuffleMove1.setConstraintMap(constraintMap);
            sim.integrators[1].getMoveManager().addMCMove(shuffleMove1);
            ((MCMoveStepTracker) shuffleMove1.getTracker()).setAcceptanceTarget(0.3);
        }

        if (false) {
            double size = (2*armLength + 5) * 1.5;
            sim.box[0].getBoundary().setBoxSize(Vector.of(size, size, size));
            sim.box[1].getBoundary().setBoxSize(Vector.of(size, size, size));
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
            DisplayBox displayBox0 = simGraphic.getDisplayBox(sim.box[0]);
            DisplayBox displayBox1 = simGraphic.getDisplayBox(sim.box[1]);
            displayBox0.setPixelUnit(new Pixel(300.0 / size));
            displayBox1.setPixelUnit(new Pixel(300.0 / size));
            displayBox0.setShowBoundary(false);
            displayBox1.setShowBoundary(false);
            ((DisplayBoxCanvasG3DSys) displayBox0.canvas).setBackgroundColor(Color.WHITE);
            ((DisplayBoxCanvasG3DSys) displayBox1.canvas).setBackgroundColor(Color.WHITE);

            ColorScheme colorScheme = new ColorScheme() {
                @Override
                public Color getAtomColor(IAtom a) {
                    Color[] c = new Color[]{Color.BLUE, Color.RED, Color.YELLOW, Color.GREEN};
                    return c[a.getParentGroup().getIndex()];
                }
            };
            displayBox0.setColorScheme(colorScheme);
            displayBox1.setColorScheme(colorScheme);

            simGraphic.makeAndDisplayFrame();

            sim.integratorOS.setNumSubSteps(1000);
            sim.setAccumulatorBlockSize(1000);

            // if running interactively, set filename to null so that it doens't read
            // (or write) to a refpref file
            ClusterAbstract myTargetCluster = targetCluster;
            Activity activityRelax = new Activity() {
                @Override
                public void runActivity(Controller.ControllerHandle handle) {
                    myTargetCluster.setTemperature(100);
                    for (int i = 0; i < 10; i++) {
                        handle.yield(sim.integrators[0]::doStep);
                        handle.yield(sim.integrators[1]::doStep);
                    }
                    myTargetCluster.setTemperature(temperature);
                }
            };
            sim.getController().addActivity(activityRelax);
            sim.initRefPref(null, 10, false);
            sim.equilibrate(null, 20, false);
            sim.getController().addActivity(new ActivityIntegrate(sim.integratorOS));
            if ((Double.isNaN(sim.refPref) || Double.isInfinite(sim.refPref) || sim.refPref == 0)) {
                throw new RuntimeException("Oops");
            }

            final DisplayTextBox stepsBox0 = new DisplayTextBox();
            stepsBox0.setLabel("Reference");
            final DisplayTextBox stepsBox1 = new DisplayTextBox();
            stepsBox1.setLabel("Target");
            JLabel jLabelPanelParentGroupSteps = new JLabel("steps");
            final JPanel panelParentGroupSteps = new JPanel(new BorderLayout());
            panelParentGroupSteps.add(jLabelPanelParentGroupSteps, CompassDirection.NORTH.toString());
            panelParentGroupSteps.add(stepsBox0.graphic(), BorderLayout.WEST);
            panelParentGroupSteps.add(stepsBox1.graphic(), BorderLayout.EAST);
            simGraphic.getPanel().controlPanel.add(panelParentGroupSteps, SimulationPanel.getVertGBC());
            DataSourceCountSteps dsSteps0 = new DataSourceCountSteps(sim.integrators[0]);
            DataPumpListener pumpSteps0 = new DataPumpListener(dsSteps0, stepsBox0, 1000);
            sim.integrators[0].getEventManager().addListener(pumpSteps0);
            DataSourceCountSteps dsSteps1 = new DataSourceCountSteps(sim.integrators[1]);
            DataPumpListener pumpSteps1 = new DataPumpListener(dsSteps1, stepsBox1, 1000);
            sim.integrators[1].getEventManager().addListener(pumpSteps1);
            final DisplayTextBox averageBox = new DisplayTextBox();
            averageBox.setLabel("Average");
            final DisplayTextBox errorBox = new DisplayTextBox();
            errorBox.setLabel("Error");
            JLabel jLabelPanelParentGroup = new JLabel("B" + nPoints);
            final JPanel panelParentGroup = new JPanel(new BorderLayout());
            panelParentGroup.add(jLabelPanelParentGroup, CompassDirection.NORTH.toString());
            panelParentGroup.add(averageBox.graphic(), BorderLayout.WEST);
            panelParentGroup.add(errorBox.graphic(), BorderLayout.EAST);
            simGraphic.getPanel().controlPanel.add(panelParentGroup, SimulationPanel.getVertGBC());

            IAction pushAnswer = new IAction() {
                DataDouble data = new DataDouble();

                public void actionPerformed() {
                    double[] ratioAndError = sim.dvo.getAverageAndError();
                    double ratio = ratioAndError[0];
                    double error = ratioAndError[1];
                    data.x = ratio;
                    averageBox.putData(data);
                    data.x = error;
                    errorBox.putData(data);
                }
            };
            IDataInfo dataInfo = new DataDouble.DataInfoDouble("B" + nPoints, new CompoundDimension(new Dimension[]{new DimensionRatio(Volume.DIMENSION, Quantity.DIMENSION)}, new double[]{nPoints - 1}));
            averageBox.putDataInfo(dataInfo);
            averageBox.setLabel("average");
            errorBox.putDataInfo(dataInfo);
            errorBox.setLabel("error");
            errorBox.setPrecision(2);
            sim.integratorOS.getEventManager().addListener(new IntegratorListenerAction(pushAnswer));
            return;
        }

        long t1 = System.nanoTime();

        targetCluster.setTemperature(temperature * 100);
        long relaxSteps = steps/40;
        Activity activityInitRefPref = new Activity() {
            @Override
            public void runActivity(Controller.ControllerHandle handle) {
                for (int i = 0; i < relaxSteps; i++) {
                    handle.yield(sim.integrators[0]::doStep);
                    handle.yield(sim.integrators[1]::doStep);
                }
            }
        };
        sim.getController().runActivityBlocking(activityInitRefPref);
        targetCluster.setTemperature(temperature);

        // if running interactively, don't use the file
        String refFileName = isCommandline ? "refpref"+nPoints+"_"+temperature : null;
        long initSteps = steps/40;
        if (params.fFile != null) {
            initSteps = steps;
        }
        // this will either read the refpref in from a file or run a short simulation to find it
        sim.initRefPref(refFileName, initSteps);

        // run another short simulation to find MC move step sizes and maybe narrow in more on the best ref pref
        // if it does continue looking for a pref, it will write the value to the file
        sim.equilibrate(refFileName, 2*initSteps);
        ActivityIntegrate ai = new ActivityIntegrate(sim.integratorOS, 1000);
        sim.setAccumulatorBlockSize(steps);

        System.out.println("equilibration finished");
        sim.setAccumulatorBlockSize(steps);
        sim.integratorOS.setNumSubSteps((int)steps);
        System.out.println("MC Move step sizes (ref)    "
                +(sim.mcMoveRotate==null ? "" : (""+sim.mcMoveRotate[0].getStepSize())));
        System.out.println("MC Move step sizes (target) "+sim.mcMoveTranslate[1].getStepSize()+" "
                +(sim.mcMoveRotate==null ? "" : (""+sim.mcMoveRotate[1].getStepSize())));

        sim.integratorOS.getMoveManager().setEquilibrating(false);

        ClusterAbstract tc = targetCluster;
        FileWriter fw;
        if (params.fFile != null) {
            try {
                fw = new FileWriter(params.fFile);
            } catch (IOException ex) {
                throw new RuntimeException(ex);
            }
            sim.integrators[1].getEventManager().addListener(new IntegratorListenerAction(new IAction() {
                @Override
                public void actionPerformed() {
                    double f = -2 * tc.value(sim.box[1]);
                    try {
                        fw.write(sim.integrators[1].getStepCount() + " " + f + "\n");
                    } catch (IOException ex) {
                        throw new RuntimeException(ex);
                    }
                }
            }));
        }
        else {
            fw = null;
        }

        sim.getController().runActivityBlocking(ai);
        long t2 = System.nanoTime();

        if (fw != null) {
            try {
                fw.close();
            } catch (IOException ex) {
                throw new RuntimeException(ex);
            }
        }

        System.out.println("final reference step frequency "+sim.integratorOS.getIdealRefStepFraction());
        System.out.println("actual reference step frequency "+sim.integratorOS.getRefStepFraction());
        sim.printResults(refIntegral);
        System.out.println("time: "+(t2-t1)/1e9);
    }

    /**
     * Inner class for parameters
     */
    public static class VirialStarParams extends ParameterBase {
        public int nPoints = 2;
        public int armLength = 3;
        public double temperature = 1;
        public long numSteps = 1000000;
        public double refFreq = -1;
        public String fFile = null;
        public int numArms = 2;
    }
}
