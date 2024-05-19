/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.config.ConformationLinear;
import etomica.graphics.*;
import etomica.potential.*;
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
import etomica.virial.MayerGeneral;
import etomica.virial.MayerTheta;
import etomica.virial.cluster.ClusterAbstract;
import etomica.virial.cluster.ClusterSum;
import etomica.virial.cluster.ClusterWeightAbs;
import etomica.virial.cluster.VirialDiagrams;
import etomica.virial.mcmove.MCMoveClusterAngle;
import etomica.virial.mcmove.MCMoveClusterMoleculeMulti;

import javax.swing.*;
import java.awt.*;
import java.util.ArrayList;
import java.util.List;

/**
 * Computes dbeta/dk for a stiff polymer chain
 */
public class VirialChainTheta {

    public static void main(String[] args) {
        VirialChainParams params = new VirialChainParams();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        } else {
            params.nPoints = 2;
            params.nSpheres = 8;
            params.temperature = 4.;
            params.numSteps = 100000000;
            params.rc = 0;
            params.kBend = 100;
        }
        final int nPoints = params.nPoints;
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

        PotentialMoleculePair pTarget = new PotentialMoleculePair(space, sm);
        System.out.println(nSpheres+"-mer chain B"+nPoints+" at T = "+temperature);
        System.out.println("Bond length: "+bondLength);
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
        MayerTheta fTheta1 = new MayerTheta(pTarget, false);
        MayerTheta fTheta2 = new MayerTheta(pTarget, true);

        boolean flex = nSpheres > 2 && nPoints > 2;
        VirialDiagrams alkaneDiagrams = new VirialDiagrams(nPoints, false, flex);
        alkaneDiagrams.setDoReeHoover(false);
        ClusterWeightAbs sampleCluster = new ClusterWeightAbs(alkaneDiagrams.makeVirialCluster(fTarget));
        ClusterSum cluster1 = alkaneDiagrams.makeVirialCluster(fTheta1);
        ClusterSum cluster2 = alkaneDiagrams.makeVirialCluster(fTheta2);
        final double refIntegral = -1;

        sampleCluster.setTemperature(temperature);
        cluster1.setTemperature(temperature);
        cluster2.setTemperature(temperature);

        // eovererr expects this string, BnHS
        System.out.println("B"+nPoints+"HS: "+refIntegral);

        PotentialMasterBonding.FullBondingInfo bondingInfo = new PotentialMasterBonding.FullBondingInfo(sm) {
            @Override
            public boolean skipBondedPair(boolean isPureAtoms, IAtom iAtom, IAtom jAtom) {
                return false;
            }
        };

        if (nSpheres < 3) {
            throw new RuntimeException("This simulation is only useful for 3+ spheres");
        }
        List<int[]> triplets = new ArrayList<>();
        for (int i=0; i<nSpheres-2; i++) {
            triplets.add(new int[]{i,i+1,i+2});
        }
        if (kBend > 0) {
            P3BondAngleStiffChain p3 = new P3BondAngleStiffChain(kBend);
            bondingInfo.setBondingPotentialTriplet(species, p3, triplets);
        }

        final SimulationVirial sim = new SimulationVirial(space, new ISpecies[]{species},
                new int[]{flex ? (nPoints+1) : nPoints}, temperature, sampleCluster, cluster2, new ClusterAbstract[]{cluster1});
        sim.setDoWiggle(false);
        sim.setBondingInfo(bondingInfo);
        sim.setIntraPairPotentials(pTarget.getAtomPotentials());
        sim.init();

        if (flex) {
            int[] constraintMap = new int[nPoints+1];
            for (int i=0; i<nPoints; i++) {
                constraintMap[i] = i;
            }
            constraintMap[nPoints] = 0;
            ((MCMoveClusterMoleculeMulti)sim.mcMoveTranslate).setConstraintMap(constraintMap);
        }

        PotentialMasterBonding.FullBondingInfo bondingInfodk = new PotentialMasterBonding.FullBondingInfo(sm);

        PotentialMasterBonding pmBondingdk = new PotentialMasterBonding(sm, sim.box, bondingInfodk);
        P3BondAngleStiffChain p3dk = new P3BondAngleStiffChain(1);
        bondingInfodk.setBondingPotentialTriplet(species, p3dk, triplets);
        fTheta1.setPotentialDK(pmBondingdk);

        fTheta2.setPotentialDK(sim.integrator.getPotentialCompute());

        System.out.println(steps+" steps");

        pTarget.setAtomPotential(type, type, p2);

        IntArrayList[] bonding = new IntArrayList[nSpheres];
        bonding[0] = new IntArrayList(new int[]{1});
        for (int i=1; i<nSpheres-1; i++) {
            bonding[i] = new IntArrayList(new int[]{i-1,i+1});
        }
        bonding[nSpheres-1] = new IntArrayList(new int[]{nSpheres-2});

        if (kBend > 0) {
            MCMoveClusterAngle angleMove = new MCMoveClusterAngle(sim.integrator.getPotentialCompute(), space, bonding, sim.getRandom(), 1);
            angleMove.setBox(sim.box);
            sim.integrator.getMoveManager().addMCMove(angleMove);
        }

        for (IAtom a : sim.box.getMoleculeList().get(1).getChildList()) {
            a.getPosition().PE(Vector.of(0,2,0));
        }
        sim.box.getMoleculeList().get(1).getChildList().get(0).getPosition().PE(Vector.of(0,0.1,0));
        sim.box.trialNotify();
        sim.box.acceptNotify();

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

        sim.equilibrate(steps/10);

        ActivityIntegrate ai = new ActivityIntegrate(sim.integrator, steps);
        sim.setAccumulatorBlockSize(steps);

        System.out.println("equilibration finished");
        sim.setAccumulatorBlockSize(steps/1000);
        System.out.println("MC Move step sizes "+sim.mcMoveTranslate.getStepSize()+" "
                +(sim.mcMoveRotate==null ? "" : (""+sim.mcMoveRotate.getStepSize())));

        sim.getController().runActivityBlocking(ai);
        long t2 = System.nanoTime();

        sim.printResults(refIntegral);
        System.out.println("time: "+(t2-t1)/1e9);
    }

    public enum TruncationChoice {
        SIMPLE, SHIFT, FORCESHIFT
    }

    /**
     * Inner class for parameters
     */
    public static class VirialChainParams extends ParameterBase {
        public int nPoints = 2;
        public int nSpheres = 3;
        public double temperature = 1;
        public long numSteps = 1000000;
        public double rc = 0;
        public double bondLength = 1;
        public TruncationChoice truncation = TruncationChoice.SHIFT;
        public double kBend = 0;
    }
}
