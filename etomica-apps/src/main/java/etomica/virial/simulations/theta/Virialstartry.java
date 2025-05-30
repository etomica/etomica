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
import etomica.potential.P2Fene;
import etomica.potential.P2WCA;
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

public class Virialstartry {

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
            params.armLength = 15;
            params.numArms = 3;
            params.temperature = 4.6;
            params.numSteps = 100000000;  // 100M steps to account for steps/=1000
        }
        final int nPoints = params.nPoints;
        int numArms = params.numArms;
        int armLength = params.armLength;
        double temperature = params.temperature;
        long steps = params.numSteps;
        double refFreq = params.refFreq;
        double sigmaHSRef = 3.5 * Math.sqrt(armLength);

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
        IPotential2 p2 = new P2WCA(1.0, 1.0);

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
        System.out.println("B"+nPoints+"HS: "+refIntegral);

        PotentialMasterBonding.FullBondingInfo bondingInfo = new PotentialMasterBonding.FullBondingInfo(sm) {
            @Override
            public boolean skipBondedPair(boolean isPureAtoms, IAtom iAtom, IAtom jAtom) {
                return false;  // Allow both FENE and WCA for bonded pairs
            }
        };

        // Use proper P2Fene class as suggested by professor
        IPotential2 pBonding = new P2Fene(1.5, 30.0);  // r0=1.5, amplitude=30.0 from paper

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

        MCMoveClusterAngle[] angleMoves = new MCMoveClusterAngle[2];
        angleMoves[0] = new MCMoveClusterAngle(pc0, space, bonding, sim.getRandom(), 1);
        angleMoves[0].setBox(sim.box[0]);
        angleMoves[0].setConstraintMap(constraintMap);
        sim.integrators[0].getMoveManager().addMCMove(angleMoves[0]);
        angleMoves[1] = new MCMoveClusterAngle(pc1, space, bonding, sim.getRandom(), 1);
        angleMoves[1].setBox(sim.box[1]);
        angleMoves[1].setConstraintMap(constraintMap);
        sim.integrators[1].getMoveManager().addMCMove(angleMoves[1]);

        // Add stretch moves for flexible FENE bonds as suggested by professor
        MCMoveClusterStretch[] stretchMoves = new MCMoveClusterStretch[2];
        stretchMoves[0] = new MCMoveClusterStretch(pc0, space, bonding, sim.getRandom(), 1);
        stretchMoves[0].setBox(sim.box[0]);
        sim.integrators[0].getMoveManager().addMCMove(stretchMoves[0]);
        stretchMoves[1] = new MCMoveClusterStretch(pc1, space, bonding, sim.getRandom(), 1);
        stretchMoves[1].setBox(sim.box[1]);
        sim.integrators[1].getMoveManager().addMCMove(stretchMoves[1]);

        if (armLength > 5) {
            MCMoveClusterShuffle shuffleMove0 = new MCMoveClusterShuffle(pc0, space, sim.getRandom());
            shuffleMove0.setBox(sim.box[0]);
            shuffleMove0.setBonding(bonding);
            shuffleMove0.setStepSizeMax(armLength - 2);
            shuffleMove0.setConstraintMap(constraintMap);
            sim.integrators[0].getMoveManager().addMCMove(shuffleMove0);
            ((MCMoveStepTracker) shuffleMove0.getTracker()).setAcceptanceTarget(0.3);
            MCMoveClusterShuffle shuffleMove1 = new MCMoveClusterShuffle(pc1, space, sim.getRandom());
            shuffleMove1.setBonding(bonding);
            shuffleMove1.setBox(sim.box[1]);
            shuffleMove1.setStepSizeMax(armLength - 2);
            shuffleMove1.setConstraintMap(constraintMap);
            sim.integrators[1].getMoveManager().addMCMove(shuffleMove1);
            ((MCMoveStepTracker) shuffleMove1.getTracker()).setAcceptanceTarget(0.3);
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

        String refFileName = isCommandline ? "refpref"+nPoints+"_"+temperature : null;
        long initSteps = steps/10;
        if (params.fFile != null) {
            initSteps = steps;
        }
        sim.initRefPref(refFileName, initSteps);

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

    public static class VirialStarParams extends ParameterBase {
        public int nPoints = 2;
        public int armLength = 15;
        public double temperature = 4.6;
        public long numSteps = 100000000;  // Increased to 100M to account for steps/=1000
        public double refFreq = -1;
        public String fFile = null;
        public int numArms = 3;
    }
}