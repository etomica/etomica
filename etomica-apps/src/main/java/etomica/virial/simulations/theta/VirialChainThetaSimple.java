/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations.theta;

import etomica.action.IAction;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.config.ConformationLinear;
import etomica.graphics.*;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveStepTracker;
import etomica.potential.*;
import etomica.potential.compute.PotentialCompute;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesBuilder;
import etomica.species.SpeciesManager;
import etomica.units.Pixel;
import etomica.util.Constants.CompassDirection;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.util.collections.IntArrayList;
import etomica.virial.BoxCluster;
import etomica.virial.MayerFunction;
import etomica.virial.MayerGeneral;
import etomica.virial.cluster.ClusterAbstract;
import etomica.virial.cluster.ClusterBonds;
import etomica.virial.cluster.ClusterSum;
import etomica.virial.cluster.ClusterWeightAbs;
import etomica.virial.mcmove.MCMoveClusterAngle;
import etomica.virial.simulations.SimulationVirial;

import javax.swing.*;
import java.awt.*;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * Computes dbeta/dk for a stiff polymer chain
 */
public class VirialChainThetaSimple {

    public static ClusterSum makeFCluster(MayerFunction f) {
        return new ClusterSum(new ClusterBonds[]{new ClusterBonds(2, new int[][][]{{{0,1}}})}, new double[]{1}, new MayerFunction[]{f});
    }

    public static void main(String[] args) {
        VirialChainParams params = new VirialChainParams();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        } else {
            params.nSpheres = 64;
            params.temperature = 4; //4.34913746887839;
            params.numSteps = 1000000;
            params.rc = 0;
            params.kBend = 0;
            params.dlnqdb = 0.7366170977353794;
            params.dlnqdk = -0.39372299306802533;
        }
        final int nPoints = 2;
        int nSpheres = params.nSpheres;
        double temperature = params.temperature;
        long steps = params.numSteps;
        double rc = params.rc;
        double bondLength = params.bondLength;
        double kBend = params.kBend;

        Space space = Space3D.getInstance();

        ConformationLinear conf = new ConformationLinear(Space3D.getInstance(), bondLength, new double[]{0,Math.PI/2});
        AtomType type = AtomType.simple("A");
        ISpecies species = new SpeciesBuilder(Space3D.getInstance())
                .addCount(type, nSpheres)
                .withConformation(conf).build();
        SpeciesManager sm = new SpeciesManager.Builder().addSpecies(species).build();

        PotentialMoleculePairCached pTarget = new PotentialMoleculePairCached(space, sm);
        System.out.println(nSpheres+"-mer chain B"+nPoints+" at T = "+temperature);
        System.out.println("Bond length: "+bondLength);
        System.out.println("kBend: "+kBend);
        if (rc>0 && rc<Double.POSITIVE_INFINITY) System.out.println("LJ truncated at "+rc+" with "+params.truncation);
        else System.out.println("LJ untruncated");
        IPotential2 p2 = new P2LennardJones(1, 1);
        if (rc > 0 && rc < Double.POSITIVE_INFINITY) {
            if (params.truncation == TruncationChoice.SIMPLE) {
                p2 = new P2SoftSphericalTruncated(p2, rc);
            }
            else if (params.truncation == TruncationChoice.SHIFT){
                p2 = new P2SoftSphericalTruncatedShifted(p2, rc);
            }
            else if (params.truncation == TruncationChoice.FORCESHIFT) {
                p2 = new P2SoftSphericalTruncatedForceShifted(p2, rc);
            }
        }

        MayerGeneral fTarget = new MayerGeneral(pTarget);
        MayerTheta11b f11b = new MayerTheta11b(pTarget);
        MayerTheta fAlpha = new MayerTheta(pTarget, true);
        MayerTheta fdfdk = new MayerTheta(pTarget, false);
        MayerTheta11a f11a = new MayerTheta11a(pTarget);
        fAlpha.setu1(-params.dlnqdb);
        fdfdk.setu1(-params.dlnqdk*temperature);

        // search needs fTarget, fAlpha, f11b
        // integration needs fAlpha, fdfdk
        // stability needs f11a, f11b, fTarget, fAlpha

        ClusterAbstract[] targetDiagrams = new ClusterAbstract[5];
        targetDiagrams[0] = makeFCluster(fTarget);
        targetDiagrams[0].setTemperature(temperature);
        targetDiagrams[1] = makeFCluster(f11b);
        targetDiagrams[1].setTemperature(temperature);
        targetDiagrams[2] = makeFCluster(fAlpha);
        targetDiagrams[2].setTemperature(temperature);
        targetDiagrams[3] = makeFCluster(fdfdk);
        targetDiagrams[3].setTemperature(temperature);
        targetDiagrams[4] = makeFCluster(f11a);
        targetDiagrams[4].setTemperature(temperature);
        ClusterWeightAbs sampleCluster = new ClusterWeightAbs(targetDiagrams[0]);
        final double refIntegral = 1;

        // eovererr expects this string, BnHS
        System.out.println("B"+nPoints+"HS: "+refIntegral);

        PotentialMasterBonding.FullBondingInfo bondingInfo = new PotentialMasterBonding.FullBondingInfo(sm) {
            @Override
            public boolean skipBondedPair(boolean isPureAtoms, IAtom iAtom, IAtom jAtom) {
                return false;
            }
        };

        if (nSpheres > 2 && kBend == 0) {
            // we need to do this to convince the system that the molecules are not rigid
            // if bondingInfo thinks molecules are rigid then intramolecular LJ will not be computed
            IPotential2 pBonding = new IPotential2() {
                @Override
                public double getRange() { return 2; }
                @Override
                public void u012add(double r2, double[] u012) { }
            };
            List<int[]> pairs = new ArrayList<>();
            for (int i=0; i<nSpheres-1; i++) {
                pairs.add(new int[]{i,i+1});
            }
            bondingInfo.setBondingPotentialPair(species, pBonding, pairs);
        }

        List<int[]> triplets = new ArrayList<>();
        for (int i = 0; i < nSpheres - 2; i++) {
            triplets.add(new int[]{i, i + 1, i + 2});
        }
        if (nSpheres > 2 && kBend < Double.POSITIVE_INFINITY) {
            P3BondAngleStiffChain p3 = new P3BondAngleStiffChain(kBend);
            bondingInfo.setBondingPotentialTriplet(species, p3, triplets);
        }
        ClusterAbstract refCluster = new ClusterAbstract() {
            @Override
            public ClusterAbstract makeCopy() {
                return null;
            }

            @Override
            public int pointCount() {
                return 2;
            }

            @Override
            public double value(BoxCluster box) {
                return 0;
            }

            @Override
            public void setTemperature(double temperature) {}
        };

        final SimulationVirial sim = new SimulationVirial(space, new ISpecies[]{species},
                new int[]{nPoints}, temperature, sampleCluster, refCluster, targetDiagrams);
        sim.setDoWiggle(false);
        sim.setBondingInfo(bondingInfo);
        sim.setIntraPairPotentials(pTarget.getAtomPotentials());
        sim.init();

        pTarget.setBox((BoxCluster) sim.box());

        PotentialMasterBonding.FullBondingInfo bondingInfodk = new PotentialMasterBonding.FullBondingInfo(sm);
        PotentialMasterBonding pmBondingdk = new PotentialMasterBonding(sm, sim.box, bondingInfodk);
        P3BondAngleStiffChain p3dk = new P3BondAngleStiffChain(1);
        bondingInfodk.setBondingPotentialTriplet(species, p3dk, triplets);

        PotentialCompute pc = sim.integrator.getPotentialCompute();

        f11b.setPotentialu(pc);
        f11a.setPotentialu(pc);
        f11a.setPotentialDK(pmBondingdk);
        fAlpha.setPotentialDK(pc);
        fdfdk.setPotentialDK(pmBondingdk);

        System.out.println(steps+" steps");

        pTarget.setAtomPotential(type, type, p2);

        IntArrayList[] bonding = new IntArrayList[nSpheres];
        bonding[0] = new IntArrayList(new int[]{1});
        for (int i=1; i<nSpheres-1; i++) {
            bonding[i] = new IntArrayList(new int[]{i-1,i+1});
        }
        bonding[nSpheres-1] = new IntArrayList(new int[]{nSpheres-2});

        MCMoveClusterAngle angleMove1 = null, angleMove2 = null, angleMove3 = null, angleMove4 = null;
        MCMoveClusterShuffle shuffleMove = null;
        MCMoveClusterReptate reptateMove = null;
        if (kBend < Double.POSITIVE_INFINITY && nSpheres >= 3) {
            IntegratorMC integrator = sim.integrator;
            if (nSpheres < 5) {
                angleMove1 = new MCMoveClusterAngle(pc, space, bonding, sim.getRandom(), 1);
                angleMove1.setBox(sim.box());
                integrator.getMoveManager().addMCMove(angleMove1);
            }
            else if (nSpheres < 9) {
                angleMove1 = new MCMoveClusterAngle(pc, space, bonding, sim.getRandom(), 1);
                angleMove1.setBox(sim.box());
                angleMove1.setAtomRange(0, nSpheres/4);
                integrator.getMoveManager().addMCMove(angleMove1, 0.5);
                angleMove2 = new MCMoveClusterAngle(pc, space, bonding, sim.getRandom(), 1);
                angleMove2.setBox(sim.box());
                angleMove2.setAtomRange(nSpheres/4, nSpheres);
                integrator.getMoveManager().addMCMove(angleMove2, 0.5);
            }
            else {
                angleMove1 = new MCMoveClusterAngle(pc, space, bonding, sim.getRandom(), 1);
                angleMove1.setBox(sim.box());
                angleMove1.setAtomRange(0, nSpheres / 8);
                integrator.getMoveManager().addMCMove(angleMove1, 0.25);
                angleMove2 = new MCMoveClusterAngle(pc, space, bonding, sim.getRandom(), 1);
                angleMove2.setBox(sim.box());
                angleMove2.setAtomRange(nSpheres / 8, nSpheres / 4);
                integrator.getMoveManager().addMCMove(angleMove2, 0.25);
                angleMove3 = new MCMoveClusterAngle(pc, space, bonding, sim.getRandom(), 1);
                angleMove3.setBox(sim.box());
                angleMove3.setAtomRange(nSpheres / 4, 3 * nSpheres / 8);
                integrator.getMoveManager().addMCMove(angleMove3, 0.25);
                angleMove4 = new MCMoveClusterAngle(pc, space, bonding, sim.getRandom(), 1);
                angleMove4.setBox(sim.box());
                angleMove4.setAtomRange(3 * nSpheres / 8, nSpheres);
                integrator.getMoveManager().addMCMove(angleMove4, 0.25);
            }
            reptateMove = new MCMoveClusterReptate(pc, space, sim.getRandom());
            reptateMove.setBox(sim.box());
            sim.integrator.getMoveManager().addMCMove(reptateMove);

            shuffleMove = new MCMoveClusterShuffle(pc, space, sim.getRandom());
            shuffleMove.setBox(sim.box());
            sim.integrator.getMoveManager().addMCMove(shuffleMove);
            ((MCMoveStepTracker)shuffleMove.getTracker()).setAcceptanceTarget(0.3);
        }

        if (false) {
            double size = (nSpheres + 5) * 1.5;
            sim.box.getBoundary().setBoxSize(Vector.of(new double[]{size, size, size}));
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
            DisplayBox displayBox0 = simGraphic.getDisplayBox(sim.box);
            displayBox0.setPixelUnit(new Pixel(300.0 / size));
            displayBox0.setShowBoundary(false);
            ((DisplayBoxCanvasG3DSys) displayBox0.canvas).setBackgroundColor(Color.WHITE);

            ColorSchemeRandomByMolecule colorScheme = new ColorSchemeRandomByMolecule(sim.getSpeciesManager(), sim.box, sim.getRandom());
            displayBox0.setColorScheme(colorScheme);

            simGraphic.makeAndDisplayFrame();

            sim.setAccumulatorBlockSize(1000);

            sim.getController().addActivity(new ActivityIntegrate(sim.integrator));

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

            return;
        }

        long t1 = System.nanoTime();

        if (params.fFile == null) {
            sim.equilibrate(steps / 10);
        }
        else {
            sim.equilibrate(steps);
        }

        ActivityIntegrate ai = new ActivityIntegrate(sim.integrator, steps);
        sim.setAccumulatorBlockSize(steps);

        System.out.println("equilibration finished");
        sim.setAccumulatorBlockSize(steps/1000);
        System.out.println("MC Move step sizes "+sim.mcMoveTranslate.getStepSize()+" "
                +(sim.mcMoveRotate==null ? "" : (""+sim.mcMoveRotate.getStepSize())));
        if (angleMove1!=null) {
            if (angleMove2 == null) {
                System.out.println("Angle move step size    " + angleMove1.getStepSize());
            } else if (angleMove3 == null) {
                System.out.println("Angle move step size    " + angleMove1.getStepSize() + " " + angleMove2.getStepSize());
            } else {
                System.out.println("Angle move step size    " + angleMove1.getStepSize() + " " + angleMove2.getStepSize() + " " + angleMove3.getStepSize() + " " + angleMove4.getStepSize());
            }
        }
        if (shuffleMove != null) System.out.println("Shuffle move step size    "+shuffleMove.getStepSize());

        ClusterAbstract tc = targetDiagrams[0];
        FileWriter fw;
        if (params.fFile != null) {
            try {
                fw = new FileWriter(params.fFile);
            } catch (IOException ex) {
                throw new RuntimeException(ex);
            }
            sim.integrator.getEventManager().addListener(new IntegratorListenerAction(new IAction() {
                @Override
                public void actionPerformed() {
                    double f = tc.value((BoxCluster)sim.box());
                    try {
                        fw.write(sim.integrator.getStepCount() + " " + f + "\n");
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

        try {
            if (fw != null) fw.close();
        }
        catch (IOException ex) {
            throw new RuntimeException(ex);
        }

        System.out.println();
        if (reptateMove!=null) System.out.println("Reptate move acceptance "+reptateMove.getTracker().acceptanceProbability());
        if (angleMove1 != null) {
            if (angleMove2 == null) {
                System.out.println("Angle move acceptance " + angleMove1.getTracker().acceptanceProbability());
            } else if (angleMove3 == null) {
                System.out.println("Angle move acceptance " + angleMove1.getTracker().acceptanceProbability() + " " + angleMove2.getTracker().acceptanceProbability());
            } else {
                System.out.println("Angle move acceptance " + angleMove1.getTracker().acceptanceProbability() + " " + angleMove2.getTracker().acceptanceProbability() + " " + angleMove3.getTracker().acceptanceProbability() + " " + angleMove4.getTracker().acceptanceProbability());
            }
            System.out.println("Shuffle move acceptance " + shuffleMove.getTracker().acceptanceProbability());
        }

        sim.printResults(refIntegral, new String[]{
                "Ibb", "alpha", "dfdk", "Ikb"
        });

        System.out.println("time: "+(t2-t1)/1e9);
    }

    public enum TruncationChoice {
        SIMPLE, SHIFT, FORCESHIFT
    }

    /**
     * Inner class for parameters
     */
    public static class VirialChainParams extends ParameterBase {
        public int nSpheres = 3;
        public double temperature = 1;
        public long numSteps = 1000000;
        public double rc = 0;
        public double bondLength = 1;
        public TruncationChoice truncation = TruncationChoice.SHIFT;
        public double kBend = 0;
        public double dlnqdb = 0;
        public double dlnqdk = 0;
        public String fFile = null;
    }
}
