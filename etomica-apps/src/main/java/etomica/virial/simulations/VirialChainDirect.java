/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.config.ConformationLinear;
import etomica.graphics.*;
import etomica.math.SpecialFunctions;
import etomica.molecule.IMoleculeList;
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
import etomica.virial.MayerFunction;
import etomica.virial.MayerGeneral;
import etomica.virial.MayerTheta;
import etomica.virial.MayerThetaWeighted;
import etomica.virial.cluster.ClusterAbstract;
import etomica.virial.cluster.ClusterBonds;
import etomica.virial.cluster.ClusterChainHS;
import etomica.virial.cluster.ClusterSum;
import etomica.virial.mcmove.MCMoveClusterAngle;
import etomica.virial.mcmove.MCMoveClusterMoleculeHSChain;

import javax.swing.*;
import java.awt.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Mayer sampling simulation for alkanes using the TraPPE-United Atoms force field.
 *   M.G. Martin and J.I. Siepmann, "Transferable Potentials for Phase
 *   Equilibria. 1. United-Atom Description of n-Alkanes," J. Phys. Chem. B
 *   102, 2569-2577 (1998)
 *   
 *   
 *   modified from VirialAlkaneFlex2 so that eovererr can be used
 *   March 2013
 */
public class VirialChainDirect {

    public static ClusterSum makeB2Cluster(MayerFunction f, double prefactor) {
        return new ClusterSum(new ClusterBonds[]{new ClusterBonds(2, new int[][][]{{{0,1}}})}, new double[]{prefactor}, new MayerFunction[]{f});
    }

    public static void main(String[] args) {
        VirialChainParams params = new VirialChainParams();
        boolean isCommandline = args.length > 0;
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        } else {
            params.nPoints = 2;
            params.nSpheres = 3;
            params.temperature = 3.244178219225989;
            params.eFENE = 0;
            params.numSteps = 1000000000L;
            params.kBend = new double[]{-1,-0.01,0,0.01,1};
            params.rc = 2;
        }
        final int nPoints = params.nPoints;
        int nSpheres = params.nSpheres;
        double temperature = params.temperature;
        long steps = params.numSteps;
        double sigmaHSRef = 1+1+params.rc;
        double eFENE = params.eFENE;
        double rc = params.rc;
        double bondLength = params.bondLength;
        double[] kBend = params.kBend;

        Space space = Space3D.getInstance();

        ConformationLinear conf = new ConformationLinear(Space3D.getInstance(), bondLength, new double[]{0,Math.PI/2});
        AtomType type = AtomType.simple("A");
        ISpecies species = new SpeciesBuilder(Space3D.getInstance())
                .addCount(type, nSpheres)
                .withConformation(conf).build();
        SpeciesManager sm = new SpeciesManager.Builder().addSpecies(species).build();

        PotentialMoleculePair pTarget = new PotentialMoleculePair(space, sm);
        System.out.println(nSpheres+"-mer chain B"+nPoints+" at T = "+temperature);
        System.out.println("Bond length: "+bondLength+"  with eFENE "+eFENE);
        System.out.println("Sampling flexible bond angles");
        System.out.println("Colecting data for bond angle bend constants: "+ Arrays.toString(kBend));
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
        MayerFunction fRefPos = new MayerFunction() {
            public void setBox(Box box) {
            }

            public double f(IMoleculeList pair, double r2, double beta) {
                return r2 < sigmaHSRef * sigmaHSRef ? 1 : 0;
            }
        };

        MayerTheta fThetadk = new MayerTheta(pTarget, false);

        int nk = 5;
        MayerThetaWeighted[] fWeighted = new MayerThetaWeighted[3*nk];
        ClusterAbstract[] targetDiagrams = new ClusterAbstract[3*nk];
        for (int i=0; i<nk; i++) {
            fWeighted[i] = new MayerThetaWeighted(fTarget);
            targetDiagrams[i] = makeB2Cluster(fWeighted[i], -0.5);
            fWeighted[nk+i] = new MayerThetaWeighted(fThetadk);
            targetDiagrams[nk + i] = makeB2Cluster(fWeighted[nk+i], 0.5);
            fWeighted[2*nk+i] = new MayerThetaWeighted(fRefPos);
            targetDiagrams[2*nk + i] = makeB2Cluster(fWeighted[2*nk+i], 1);
        }
        ClusterChainHS refCluster = new ClusterChainHS(nPoints, fRefPos);
        double vhs = (4.0 / 3.0) * Math.PI * sigmaHSRef * sigmaHSRef * sigmaHSRef;
        final double refIntegral = SpecialFunctions.factorial(nPoints) / 2 * Math.pow(vhs, nPoints - 1);

        refCluster.setTemperature(temperature);
        for (int i=0; i<targetDiagrams.length; i++) {
            targetDiagrams[i].setTemperature(temperature);
        }

        System.out.println("sigmaHSRef: "+sigmaHSRef);
        // eovererr expects this string, BnHS
        System.out.println("B"+nPoints+"HS: "+refIntegral);

        PotentialMasterBonding.FullBondingInfo bondingInfo = new PotentialMasterBonding.FullBondingInfo(sm) {
            @Override
            public boolean skipBondedPair(boolean isPureAtoms, IAtom iAtom, IAtom jAtom) {
                return false;
            }
        };

        IPotential2 pBonding;
        if (eFENE > 0 && eFENE < Double.POSITIVE_INFINITY) {
            pBonding = new P2Fene(2, eFENE);
        }
        else {
            // we need to do this to convince the system that the molecules are not rigid
            // if bondingInfo thinks molecules are rigid then intramolecular LJ will not be computed
            pBonding = new IPotential2() {
                @Override
                public double getRange() { return 2; }
                @Override
                public void u012add(double r2, double[] u012) { }
            };
        }
        List<int[]> pairs = new ArrayList<>();
        for (int i=0; i<nSpheres-1; i++) {
            pairs.add(new int[]{i,i+1});
        }
        bondingInfo.setBondingPotentialPair(species, pBonding, pairs);
        List<int[]> triplets = new ArrayList<>();
        for (int i=0; i<nSpheres-2; i++) {
            triplets.add(new int[]{i,i+1,i+2});
        }

        P3BondAngleStiffChain p3 = new P3BondAngleStiffChain(0);
        bondingInfo.setBondingPotentialTriplet(species, p3, triplets);

        final SimulationVirial sim = new SimulationVirial(space, new ISpecies[]{species}, new int[]{nPoints},
                temperature, refCluster, refCluster, targetDiagrams);
        sim.setDoWiggle(false);
        sim.setBondingInfo(bondingInfo);
        sim.setIntraPairPotentials(pTarget.getAtomPotentials());
        sim.init();

        PotentialMasterBonding.FullBondingInfo bondingInfodk0 = new PotentialMasterBonding.FullBondingInfo(sm);

        PotentialMasterBonding pmBondingdk0 = new PotentialMasterBonding(sm, sim.box, bondingInfodk0);
        P3BondAngleStiffChain p3dk0 = new P3BondAngleStiffChain(1);
        bondingInfodk0.setBondingPotentialTriplet(species, p3dk0, triplets);
        fThetadk.setPotentialDK(pmBondingdk0);
        for (int i=0; i<nk; i++) {
            PotentialMasterBonding.FullBondingInfo bondingInfodk = new PotentialMasterBonding.FullBondingInfo(sm);
            PotentialMasterBonding pmBondingdk = new PotentialMasterBonding(sm, sim.box, bondingInfodk);
            P3BondAngleStiffChain p3dk = new P3BondAngleStiffChain(kBend[i]);
            bondingInfodk.setBondingPotentialTriplet(species, p3dk, triplets);
            fWeighted[i].setPotentialExtra(pmBondingdk);
            fWeighted[nk+i].setPotentialExtra(pmBondingdk);
            fWeighted[2*nk+i].setPotentialExtra(pmBondingdk);
        }

        sim.integrator.getMoveManager().removeMCMove(sim.mcMoveTranslate);
        MCMoveClusterMoleculeHSChain mcMoveHSC = new MCMoveClusterMoleculeHSChain(sim.getRandom(), sim.box, sigmaHSRef);
        sim.integrator.getMoveManager().addMCMove(mcMoveHSC);

        System.out.println(steps+" steps");

        pTarget.setAtomPotential(type, type, p2);

        IntArrayList[] bonding = new IntArrayList[nSpheres];
        bonding[0] = new IntArrayList(new int[]{1});
        for (int i=1; i<nSpheres-1; i++) {
            bonding[i] = new IntArrayList(new int[]{i-1,i+1});
        }
        bonding[nSpheres-1] = new IntArrayList(new int[]{nSpheres-2});

        MCMoveClusterAngle angleMove = new MCMoveClusterAngle(sim.integrator.getPotentialCompute(), space, bonding, sim.getRandom(), 1);
        angleMove.setBox(sim.box);
        sim.integrator.getMoveManager().addMCMove(angleMove);

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

            // if running interactively, set filename to null so that it doens't read
            // (or write) to a refpref file
//            sim.equilibrate(1000000);
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
        sim.setAccumulatorBlockSize(steps/1000);

        System.out.println("equilibration finished");
        System.out.println("MC Move step sizes (ref)    "
                +(sim.mcMoveRotate==null ? "" : (""+sim.mcMoveRotate.getStepSize())));

        sim.integrator.getMoveManager().setEquilibrating(false);
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
        public double eFENE = 0;
        public double rc = 2;
        public double bondLength = 1;
        public TruncationChoice truncation = TruncationChoice.SIMPLE;
        public double[] kBend = new double[]{0};
    }
}
