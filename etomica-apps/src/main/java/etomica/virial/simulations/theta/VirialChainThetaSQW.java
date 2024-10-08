/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations.theta;

import etomica.action.IAction;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.DiameterHashByType;
import etomica.config.ConformationLinear;
import etomica.graphics.*;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveStepTracker;
import etomica.potential.IPotential2;
import etomica.potential.P2HardGeneric;
import etomica.potential.P2SquareWell;
import etomica.potential.PotentialMasterBonding;
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
import etomica.virial.mcmove.MCMoveClusterMoleculeMulti;
import etomica.virial.mcmove.MCMoveClusterRotateMoleculeMulti;
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
public class VirialChainThetaSQW {

    public static ClusterSum makeFCluster(MayerFunction f) {
        return new ClusterSum(new ClusterBonds[]{new ClusterBonds(2, new int[][][]{{{0,1}}})}, new double[]{1}, new MayerFunction[]{f});
    }

    public static void main(String[] args) {
        VirialChainParams params = new VirialChainParams();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        } else {
            params.nH = 8;
            params.temperature = 5.8;
            params.numSteps = 10000000;
            params.lattice = true;
        }
        final int nPoints = 2;
        int nH = params.nH;
        int nT = params.nT;
        double lambda = params.lambda;
        double temperature = params.temperature;
        long steps = params.numSteps;
        int nSpheres = nH + nT;
        boolean doLattice = params.lattice;
        int latNabs = params.latNabs;

        double epsH = 1;
        double epsT = 3;
        double epsHT = 1;
        System.out.println("epsH: "+epsH+"   epsT: "+epsT+"    epsHT: "+epsHT);

        double sigma = doLattice ? 0.9 : 1;
        if (doLattice) {
            if (latNabs == 26) lambda = 2;
            else if (latNabs == 20) lambda = 1.6;
            else if (latNabs == 6) lambda = 1.2;
        }

        Space space = Space3D.getInstance();

        ConformationLinear conf = new ConformationLinear(Space3D.getInstance(), 1, new double[]{0,Math.PI/2});
        AtomType typeH = AtomType.simple("H");
        AtomType typeT = AtomType.simple("T");
        ISpecies species = new SpeciesBuilder(Space3D.getInstance())
                .addCount(typeH, nH)
                .addCount(typeT, nT)
                .withConformation(conf).build();
        SpeciesManager sm = new SpeciesManager.Builder().addSpecies(species).build();

        PotentialMoleculePairCached pTarget = new PotentialMoleculePairCached(space, sm);
        System.out.println(nSpheres+"-mer chain");
        System.out.println("T: "+temperature);
        System.out.println("nH: "+nH);
        System.out.println("nT: "+nT);
        if (!doLattice) System.out.println("SQW lambda: "+lambda);
        P2HardGeneric p2HH = P2SquareWell.makePotential(sigma, lambda, epsH);
        P2HardGeneric p2TT = P2SquareWell.makePotential(sigma, lambda, epsT);
        P2HardGeneric p2HT = P2SquareWell.makePotential(sigma, lambda, epsHT);
        if (nH>0) pTarget.setAtomPotential(typeH, typeH, p2HH);
        if (nT>0) pTarget.setAtomPotential(typeT, typeT, p2TT);
        if (nH*nT>0) pTarget.setAtomPotential(typeH, typeT, p2HT);

        MayerGeneral fTarget = new MayerGeneral(pTarget);
        MayerTheta11b f11b = new MayerTheta11b(pTarget);
        MayerTheta fAlpha = new MayerTheta(pTarget, true);
        fAlpha.setu1(-params.dlnqdb);

        // search needs fTarget, fAlpha, f11b

        ClusterAbstract[] targetDiagrams = new ClusterAbstract[3];
        targetDiagrams[0] = makeFCluster(fTarget);
        targetDiagrams[0].setTemperature(temperature);
        targetDiagrams[1] = makeFCluster(f11b);
        targetDiagrams[1].setTemperature(temperature);
        targetDiagrams[2] = makeFCluster(fAlpha);
        targetDiagrams[2].setTemperature(temperature);
        ClusterWeightAbs sampleCluster = new ClusterWeightAbs(targetDiagrams[0]);
        final double refIntegral = 1;

        // eovererr expects this string, BnHS
        System.out.println("B"+nPoints+"HS: "+refIntegral);

        PotentialMasterBonding.FullBondingInfo bondingInfo = new PotentialMasterBonding.FullBondingInfo(sm);

        if (nSpheres > 2) {
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

        PotentialCompute pc = sim.integrator.getPotentialCompute();

        f11b.setPotentialu(pc);
        fAlpha.setPotentialDK(pc);

        System.out.println(steps+" steps");

        IntArrayList[] bonding = new IntArrayList[nSpheres];
        bonding[0] = new IntArrayList(new int[]{1});
        for (int i=1; i<nSpheres-1; i++) {
            bonding[i] = new IntArrayList(new int[]{i-1,i+1});
        }
        bonding[nSpheres-1] = new IntArrayList(new int[]{nSpheres-2});

        if (nSpheres==1) {
            sim.integrator.getMoveManager().removeMCMove(sim.mcMoveTranslate);
            sim.mcMoveTranslate = new MCMoveClusterMoleculeMulti(sim.getRandom(), sim.box());
            sim.integrator.getMoveManager().addMCMove(sim.mcMoveTranslate);
        }
        ((MCMoveClusterMoleculeMulti)sim.mcMoveTranslate).setDoLattice(doLattice);
        ((MCMoveClusterRotateMoleculeMulti)sim.mcMoveRotate).setDoLattice(doLattice);

        MCMoveClusterAngle angleMove1 = null, angleMove2 = null, angleMove3 = null, angleMove4 = null;
        MCMoveClusterShuffle shuffleMove = null;
        MCMoveClusterReptate reptateMove = null;
        if (nSpheres >= 3) {
            IntegratorMC integrator = sim.integrator;
            if (nSpheres < 6 || doLattice) {
                angleMove1 = new MCMoveClusterAngle(pc, space, bonding, sim.getRandom(), 1);
                angleMove1.setBox(sim.box());
                angleMove1.setDoLattice(doLattice);
                integrator.getMoveManager().addMCMove(angleMove1);
            }
            else if (nSpheres < 9) {
                angleMove1 = new MCMoveClusterAngle(pc, space, bonding, sim.getRandom(), 1);
                angleMove1.setBox(sim.box());
                angleMove1.setAtomRange(0, nSpheres/4);
                angleMove1.setDoLattice(doLattice);
                integrator.getMoveManager().addMCMove(angleMove1, 0.5);
                angleMove2 = new MCMoveClusterAngle(pc, space, bonding, sim.getRandom(), 1);
                angleMove2.setBox(sim.box());
                angleMove2.setAtomRange(nSpheres/4, nSpheres);
                angleMove2.setDoLattice(doLattice);
                integrator.getMoveManager().addMCMove(angleMove2, 0.5);
            }
            else {
                angleMove1 = new MCMoveClusterAngle(pc, space, bonding, sim.getRandom(), 1);
                angleMove1.setBox(sim.box());
                angleMove1.setAtomRange(0, nSpheres / 8);
                angleMove1.setDoLattice(doLattice);
                integrator.getMoveManager().addMCMove(angleMove1, 0.25);
                angleMove2 = new MCMoveClusterAngle(pc, space, bonding, sim.getRandom(), 1);
                angleMove2.setBox(sim.box());
                angleMove2.setAtomRange(nSpheres / 8, nSpheres / 4);
                angleMove2.setDoLattice(doLattice);
                integrator.getMoveManager().addMCMove(angleMove2, 0.25);
                angleMove3 = new MCMoveClusterAngle(pc, space, bonding, sim.getRandom(), 1);
                angleMove3.setBox(sim.box());
                angleMove3.setAtomRange(nSpheres / 4, 3 * nSpheres / 8);
                angleMove3.setDoLattice(doLattice);
                integrator.getMoveManager().addMCMove(angleMove3, 0.25);
                angleMove4 = new MCMoveClusterAngle(pc, space, bonding, sim.getRandom(), 1);
                angleMove4.setBox(sim.box());
                angleMove4.setAtomRange(3 * nSpheres / 8, nSpheres);
                angleMove4.setDoLattice(doLattice);
                integrator.getMoveManager().addMCMove(angleMove4, 0.25);
            }
            reptateMove = new MCMoveClusterReptate(pc, space, sim.getRandom());
            reptateMove.setBox(sim.box());
            if (doLattice && latNabs == 6) {
                reptateMove.setDoCubicLattice(true);
            }
            else {
                reptateMove.setDoLattice(doLattice);
            }
            sim.integrator.getMoveManager().addMCMove(reptateMove);

            if (nSpheres > 3) {
                shuffleMove = new MCMoveClusterShuffle(pc, space, sim.getRandom());
                shuffleMove.setBox(sim.box());
                shuffleMove.setDoLattice(doLattice);
                sim.integrator.getMoveManager().addMCMove(shuffleMove);
                ((MCMoveStepTracker) shuffleMove.getTracker()).setAcceptanceTarget(0.3);
            }
        }

        if (false) {
            double size = (nSpheres + 5) * 1.5;
            sim.box.getBoundary().setBoxSize(Vector.of(new double[]{size, size, size}));
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
            DisplayBox displayBox0 = simGraphic.getDisplayBox(sim.box);
            displayBox0.setPixelUnit(new Pixel(300.0 / size));
            displayBox0.setShowBoundary(false);
            ((DisplayBoxCanvasG3DSys) displayBox0.canvas).setBackgroundColor(Color.WHITE);
            ((DiameterHashByType)displayBox0.getDiameterHash()).setDiameter(typeH, 1);
            ((DiameterHashByType)displayBox0.getDiameterHash()).setDiameter(typeT, 1);

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
        if (doLattice) {
            System.out.println("MC Move step size " + sim.mcMoveTranslate.getStepSize());
        }
        else {
            System.out.println("MC Move step sizes " + sim.mcMoveTranslate.getStepSize() + " "
                    + (sim.mcMoveRotate == null ? "" : ("" + sim.mcMoveRotate.getStepSize())));
        }
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
        if (doLattice) {
            System.out.println("Translate move acceptance "+sim.mcMoveTranslate.getTracker().acceptanceProbability());
            System.out.println("Rotate move acceptance "+sim.mcMoveRotate.getTracker().acceptanceProbability());
        }
        if (reptateMove!=null) System.out.println("Reptate move acceptance "+reptateMove.getTracker().acceptanceProbability());
        if (angleMove1 != null) {
            if (angleMove2 == null) {
                System.out.println("Angle move acceptance " + angleMove1.getTracker().acceptanceProbability());
            } else if (angleMove3 == null) {
                System.out.println("Angle move acceptance " + angleMove1.getTracker().acceptanceProbability() + " " + angleMove2.getTracker().acceptanceProbability());
            } else {
                System.out.println("Angle move acceptance " + angleMove1.getTracker().acceptanceProbability() + " " + angleMove2.getTracker().acceptanceProbability() + " " + angleMove3.getTracker().acceptanceProbability() + " " + angleMove4.getTracker().acceptanceProbability());
            }
            if (shuffleMove!=null) System.out.println("Shuffle move acceptance " + shuffleMove.getTracker().acceptanceProbability());
        }

        sim.printResults(refIntegral, new String[]{
                "Ibb", "alpha"
        });

        System.out.println("time: "+(t2-t1)/1e9);
    }

    /**
     * Inner class for parameters
     */
    public static class VirialChainParams extends ParameterBase {
        public int nH = 5;
        public int nT = 0;
        public double temperature = 1;
        public double lambda = 1.5;
        public long numSteps = 1000000;
        public double dlnqdb = 0;
        public String fFile = null;
        public boolean lattice = false;
        public int latNabs = 26;
    }
}
