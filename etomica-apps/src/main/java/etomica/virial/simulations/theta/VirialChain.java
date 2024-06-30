/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations.theta;

import etomica.action.IAction;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.config.ConformationLinear;
import etomica.data.IDataInfo;
import etomica.data.types.DataDouble;
import etomica.graphics.*;
import etomica.integrator.IntegratorListenerAction;
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
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.*;
import etomica.util.Constants.CompassDirection;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.util.collections.IntArrayList;
import etomica.virial.MayerFunction;
import etomica.virial.MayerGeneral;
import etomica.virial.cluster.*;
import etomica.virial.mcmove.MCMoveClusterAngle;
import etomica.virial.mcmove.MCMoveClusterMoleculeHSChain;
import etomica.virial.mcmove.MCMoveClusterStretch;
import etomica.virial.simulations.SimulationVirialOverlap2;
import etomica.virial.wheatley.ClusterWheatleySoftDerivatives;

import javax.swing.*;
import java.awt.*;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
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
public class VirialChain {

    public static ClusterSum makeB2Cluster(MayerFunction f) {
        return new ClusterSum(new ClusterBonds[]{new ClusterBonds(2, new int[][][]{{{0,1}}})}, new double[]{1}, new MayerFunction[]{f});
    }

    public static void main(String[] args) {
        VirialChainParams params = new VirialChainParams();
        boolean isCommandline = args.length > 0;
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        } else {
            params.nPoints = 4;
            params.nSpheres = 128;
            params.temperature = 4.55897667284274;
            params.kBend = 7.11111111111111111;
            params.numSteps = 1000000;
        }
        final int nPoints = params.nPoints;
        int nSpheres = params.nSpheres;
        double temperature = params.temperature;
        long steps = params.numSteps;
        double refFreq = params.refFreq;
        double sigmaHSRef = 1.5 + 0.15*nSpheres;
        double eFENE = params.eFENE;
        double rc = params.rc;
        double bondLength = params.bondLength;
        double kBend = params.kBend;
        int nDer = params.nDer;
        TargetChoice targetChoice = params.targetChoice;
        if (targetChoice != TargetChoice.NORMAL && nDer > 0) {
            throw new RuntimeException("Arbitrary derivatives only available for NORMAL");
        }

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
        if (kBend == Double.POSITIVE_INFINITY) System.out.println("Rigid bond angles");
        else if (kBend == 0) System.out.println("Flexible bond angles");
        else System.out.println("Bond angle bend constant: "+kBend);
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

        boolean flex = nSpheres > 2 && nPoints > 2;
        VirialDiagrams alkaneDiagrams = new VirialDiagrams(nPoints, false, flex);
        alkaneDiagrams.setDoReeHoover(false);
        ClusterAbstract targetCluster = null;
        ClusterAbstract[] targetDiagrams = new ClusterAbstract[0];
        MayerTheta fThetadk = new MayerTheta(pTarget, false);
        MayerTheta fThetadBeta = new MayerTheta(pTarget, true);

        if (targetChoice == TargetChoice.NORMAL) {
            if (nDer == 0) {
                targetCluster = alkaneDiagrams.makeVirialCluster(fTarget);
            } else if (flex) {
                throw new RuntimeException("Cannot compute derivatives for flex diagrams");
            } else {
                targetCluster = new ClusterWheatleySoftDerivatives(nPoints, fTarget, 1e-12, nDer);
                targetDiagrams = new ClusterAbstract[nDer];
                for(int m=1;m<=nDer;m++){
                    targetDiagrams[m-1]= new ClusterWheatleySoftDerivatives.ClusterRetrievePrimes((ClusterWheatleySoftDerivatives) targetCluster,m);
                }
            }
        }
        else if (targetChoice == TargetChoice.DBETADK) {
            targetCluster = new ClusterWheatleySoftDerivatives(nPoints, fTarget, 1e-12, 1);
            targetDiagrams = new ClusterAbstract[]{makeB2Cluster(fThetadk),makeB2Cluster(fThetadBeta),
                    new ClusterWheatleySoftDerivatives.ClusterRetrievePrimes((ClusterWheatleySoftDerivatives) targetCluster,1)};
        }
        ClusterChainHS refCluster = new ClusterChainHS(nPoints, fRefPos);
        double vhs = (4.0 / 3.0) * Math.PI * sigmaHSRef * sigmaHSRef * sigmaHSRef;
        final double refIntegral = SpecialFunctions.factorial(nPoints) / 2 * Math.pow(vhs, nPoints - 1);

        targetCluster.setTemperature(temperature);
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

        if (nSpheres > 1 && (kBend == 0 || eFENE > 0)) {
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
        }
        List<int[]> triplets = new ArrayList<>();
        for (int i=0; i<nSpheres-2; i++) {
            triplets.add(new int[]{i,i+1,i+2});
        }
        if (nSpheres > 2 && kBend < Double.POSITIVE_INFINITY) {
            P3BondAngleStiffChain p3 = new P3BondAngleStiffChain(kBend);
            bondingInfo.setBondingPotentialTriplet(species, p3, triplets);
        }

        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space, sm,
                new int[]{flex ? (nPoints+1) : nPoints},temperature, refCluster, targetCluster);
        sim.setExtraTargetClusters(targetDiagrams);
        sim.setDoWiggle(false);
        sim.setBondingInfo(bondingInfo);
        sim.setIntraPairPotentials(pTarget.getAtomPotentials());
        sim.init();

        PotentialMasterBonding.FullBondingInfo bondingInfodk = new PotentialMasterBonding.FullBondingInfo(sm);

        PotentialMasterBonding pmBondingdk = new PotentialMasterBonding(sm, sim.box[1], bondingInfodk);
        P3BondAngleStiffChain p3dk = new P3BondAngleStiffChain(1);
        bondingInfodk.setBondingPotentialTriplet(species, p3dk, triplets);
        fThetadk.setPotentialDK(pmBondingdk);

        fThetadBeta.setPotentialDK(sim.integrators[1].getPotentialCompute());

        sim.integrators[0].getMoveManager().removeMCMove(sim.mcMoveTranslate[0]);
        MCMoveClusterMoleculeHSChain mcMoveHSC = new MCMoveClusterMoleculeHSChain(sim.getRandom(), sim.box[0], sigmaHSRef);
        sim.integrators[0].getMoveManager().addMCMove(mcMoveHSC);
        sim.accumulators[0].setBlockSize(1);

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

        IntArrayList[] bonding = new IntArrayList[nSpheres];
        bonding[0] = new IntArrayList(new int[]{1});
        for (int i=1; i<nSpheres-1; i++) {
            bonding[i] = new IntArrayList(new int[]{i-1,i+1});
        }
        bonding[nSpheres-1] = new IntArrayList(new int[]{nSpheres-2});

        if (nSpheres > 1 && eFENE > 0 && eFENE < Double.POSITIVE_INFINITY) {
            stretchMoves = new MCMoveClusterStretch[2];
            stretchMoves[0] = new MCMoveClusterStretch(sim.integrators[0].getPotentialCompute(), space, bonding, sim.getRandom(), 1);
            stretchMoves[0].setBox(sim.box[0]);
            sim.integrators[0].getMoveManager().addMCMove(stretchMoves[0]);
            stretchMoves[1] = new MCMoveClusterStretch(sim.integrators[1].getPotentialCompute(), space, bonding, sim.getRandom(), 1);
            stretchMoves[1].setBox(sim.box[1]);
            sim.integrators[1].getMoveManager().addMCMove(stretchMoves[1]);
        }

        if (nSpheres > 2 && kBend < Double.POSITIVE_INFINITY) {
            angleMoves = new MCMoveClusterAngle[2];
            angleMoves[0] = new MCMoveClusterAngle(sim.integrators[0].getPotentialCompute(), space, bonding, sim.getRandom(), 1);
            angleMoves[0].setBox(sim.box[0]);
            sim.integrators[0].getMoveManager().addMCMove(angleMoves[0]);
//            ((MCMoveStepTracker)angleMoves[0].getTracker()).setNoisyAdjustment(true);
            angleMoves[1] = new MCMoveClusterAngle(sim.integrators[1].getPotentialCompute(), space, bonding, sim.getRandom(), 1);
            angleMoves[1].setBox(sim.box[1]);
            sim.integrators[1].getMoveManager().addMCMove(angleMoves[1]);
//            ((MCMoveStepTracker)angleMoves[1].getTracker()).setNoisyAdjustment(true);
        }

        if(false) {
            double size = (nSpheres + 5) * 1.5;
            sim.box[0].getBoundary().setBoxSize(Vector.of(new double[]{size, size, size}));
            sim.box[1].getBoundary().setBoxSize(Vector.of(new double[]{size, size, size}));
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
            DisplayBox displayBox0 = simGraphic.getDisplayBox(sim.box[0]);
            DisplayBox displayBox1 = simGraphic.getDisplayBox(sim.box[1]);
            displayBox0.setPixelUnit(new Pixel(300.0 / size));
            displayBox1.setPixelUnit(new Pixel(300.0 / size));
            displayBox0.setShowBoundary(false);
            displayBox1.setShowBoundary(false);
            ((DisplayBoxCanvasG3DSys) displayBox0.canvas).setBackgroundColor(Color.WHITE);
            ((DisplayBoxCanvasG3DSys) displayBox1.canvas).setBackgroundColor(Color.WHITE);


            ColorSchemeRandomByMolecule colorScheme = new ColorSchemeRandomByMolecule(sim.getSpeciesManager(), sim.box[0], sim.getRandom());
            displayBox0.setColorScheme(colorScheme);
            colorScheme = new ColorSchemeRandomByMolecule(sim.getSpeciesManager(), sim.box[1], sim.getRandom());
            displayBox1.setColorScheme(colorScheme);


            simGraphic.makeAndDisplayFrame();

            sim.integratorOS.setNumSubSteps(1000);
            sim.setAccumulatorBlockSize(1000);

            // if running interactively, set filename to null so that it doens't read
            // (or write) to a refpref file
            sim.initRefPref(null, 10, false);
            sim.equilibrate(null, 20, false);
            sim.getController().addActivity(new ActivityIntegrate(sim.integratorOS));
            if ((Double.isNaN(sim.refPref) || Double.isInfinite(sim.refPref) || sim.refPref == 0)) {
                throw new RuntimeException("Oops");
            }

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
        String[] extraNames = new String[targetDiagrams.length];
        if (nDer > 0) {
            for (int i = 1; i <= nDer; i++) {
                extraNames[i - 1] = "derivative " + i;
            }
        }
        else if (targetChoice == TargetChoice.DBETADK) {
            extraNames[0] = "dB2dk";
            extraNames[1] = "dB2dbeta";
            extraNames[2] = "CWSD dB2dbeta";
        }
        sim.printResults(refIntegral, extraNames);
        System.out.println("time: "+(t2-t1)/1e9);
    }

    public enum TruncationChoice {
        SIMPLE, SHIFT, FORCESHIFT
    }

    public enum TargetChoice {
        NORMAL, DBETADK
    }

    /**
     * Inner class for parameters
     */
    public static class VirialChainParams extends ParameterBase {
        public int nPoints = 2;
        public int nSpheres = 3;
        public double temperature = 1;
        public long numSteps = 1000000;
        public double refFreq = -1;
        public double eFENE = 0;
        public double rc = 0;
        public double bondLength = 1;
        public TruncationChoice truncation = TruncationChoice.SHIFT;
        public double kBend = 0;
        public int nDer = 0;
        public TargetChoice targetChoice = TargetChoice.NORMAL;
        public String fFile = null;
    }
}
