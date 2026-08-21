/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations.theta;

// CHANGE 1: removed XYZWriter import (no xyz output needed here)
// CHANGE 2: added MeterArmProfile import
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import java.io.FileWriter;
import etomica.atom.DiameterHashByType;
import etomica.atom.IAtom;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.DataPumpListener;
import etomica.data.DataSourceCountSteps;
import etomica.data.IData;
import etomica.graphics.*;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveStepTracker;
import etomica.potential.IPotential2;
import etomica.potential.P2LennardJones;
import etomica.potential.PotentialMasterBonding;
import etomica.potential.PotentialMoleculePair;
import etomica.potential.compute.*;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularNonperiodic;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesBuilder;
import etomica.species.SpeciesManager;
import etomica.starpolymer.ConformationStarPolymerGraft;
import etomica.units.Pixel;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.util.collections.IntArrayList;
import etomica.virial.mcmove.MCMoveClusterAngle;
import etomica.virial.mcmove.MCMoveClusterAngleMulti;
import etomica.virial.mcmove.MCMoveClusterStretch;
import etomica.virial.simulations.theta.MCMoveClusterShuffle;

import java.awt.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

// CHANGE 3: class renamed from VirialStarSingle to VirialStarProfile
public class VirialStarProfile {

    public static void main(String[] args) {
        VirialStarParams params = new VirialStarParams();

        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        } else {
            params.armLength =2 ;
            params.numArms = 1;
            params.temperature = 1.5;
            params.numSteps = 10000;
            params.coreSigma = 16.0;
            params.epsCore = 1.0;
            params.ideal = false;
        }
        int numArms = params.numArms;
        int armLength = params.armLength;
        double temperature = params.temperature;
        long steps = params.numSteps;
        boolean ideal = params.ideal;
        if (armLength < 2){
            throw new RuntimeException("arm length must be at least 2");
        }
        Space space = Space3D.getInstance();

        ConformationStarPolymerGraft conf = new ConformationStarPolymerGraft(Space3D.getInstance(), numArms, armLength);
        conf.setSigma0(params.coreSigma);
        AtomType typeCore = AtomType.simple("Core");
        AtomType typeMono = AtomType.simple("Mono");
        ISpecies species = new SpeciesBuilder(Space3D.getInstance())
                .addCount(typeCore, 1)
                .addCount(typeMono, numArms*armLength)
                .withConformation(conf)
                .build();
        SpeciesManager sm = new SpeciesManager.Builder().addSpecies(species).build();

        PotentialMoleculePair pTarget = new PotentialMoleculePair(space, sm);
        System.out.println(numArms+" arms of length "+armLength+" at T = "+temperature);
        double sigC = params.coreSigma;
        double epsC = params.epsCore;
        double sigM = 1;
        double epsM = 1;
        double sigCM = 0.5*(sigC + sigM);
        double epsCM = Math.sqrt(epsC * epsM);

        IPotential2 p2CoreCore = new P2LennardJones(sigC,  epsC);
        IPotential2 p2MonoMono = new P2LennardJones(sigM,  epsM);
        IPotential2 p2CoreMono = new P2LennardJones(sigCM, epsCM);

        PotentialMasterBonding.FullBondingInfo bondingInfo = new PotentialMasterBonding.FullBondingInfo(sm) {
            @Override
            public boolean skipBondedPair(boolean isPureAtoms, IAtom iAtom, IAtom jAtom) {
                return false;
            }
        };

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

        Simulation sim = new Simulation(Space3D.getInstance(), sm);
        sim.makeBox(new BoundaryRectangularNonperiodic(space));
        sim.box().setNMolecules(species, 1);

        System.out.println("random seeds: "+ Arrays.toString(sim.getRandomSeeds()));

        pTarget.setAtomPotential(typeCore, typeCore, p2CoreCore);
        pTarget.setAtomPotential(typeMono, typeMono, p2MonoMono);
        pTarget.setAtomPotential(typeCore, typeMono, p2CoreMono);
        pTarget.setAtomPotential(typeMono, typeCore, p2CoreMono);

        PotentialMasterBonding pmBonding = new PotentialMasterBonding(sm, sim.box(), bondingInfo);
        NeighborManager nbrManager = new NeighborManagerIntra(sim.box(), bondingInfo);
        if (ideal) {
            nbrManager = new NeighborManagerIntra(sim.box(), bondingInfo) {
                public NeighborIterator makeNeighborIterator() {
                    return new NeighborIterator() {
                        @Override
                        public void iterUpNeighbors(int i, NeighborConsumer consumer) {}
                        @Override
                        public void iterDownNeighbors(int i, NeighborConsumer consumer) {}
                        @Override
                        public void iterAllNeighbors(int i, NeighborConsumer consumer) {}
                    };
                }
            };
        }
        PotentialComputePair pcPair = new PotentialComputePair(sm, sim.box(), nbrManager, pTarget.getAtomPotentials());
        PotentialComputeAggregate pc = new PotentialComputeAggregate(pmBonding, pcPair);

        IntegratorMC integrator = new IntegratorMC(pc, sim.getRandom(), temperature, sim.box());

        System.out.println(steps+" steps");

        IntArrayList[] bonding = new IntArrayList[1+numArms*armLength];
        bonding[0] = new IntArrayList(numArms);

        for (int i=0; i<numArms; i++) {
            bonding[0].add(1+i*armLength);
            if (armLength > 1) {
                bonding[1 + i * armLength] = new IntArrayList(new int[]{0, 1 + i * armLength + 1});
            } else {
                bonding[1 + i * armLength] = new IntArrayList(new int[]{0});
            }
            for (int j = 1; j < armLength-1; j++) {
                bonding[1+armLength*i+j] = new IntArrayList(new int[]{1+armLength*i+j-1, 1+armLength*i+j+1});
            }
            if (armLength > 1) {
                bonding[1 + (i+1) * armLength - 1] = new IntArrayList(new int[]{1 + (i+1) * armLength - 2});
            }
        }

        MCMoveClusterAngle angleMoveInner=null;
        if (armLength>2){
            angleMoveInner = new MCMoveClusterAngle(pc, space, bonding, sim.getRandom(), 1,
                    new MCMoveClusterAngle.AtomChooserStarfl(sim.getRandom(), numArms, armLength, 2, (armLength+1)/2));
            angleMoveInner.setBox(sim.box());
            integrator.getMoveManager().addMCMove(angleMoveInner, 0.1);
        }
        MCMoveClusterAngle angleMoveOuter = new MCMoveClusterAngle(pc, space, bonding, sim.getRandom(), 1,
                new MCMoveClusterAngle.AtomChooserStarfl(sim.getRandom(), numArms, armLength, (armLength+1) / 2 + 1, armLength));
        angleMoveOuter.setBox(sim.box());
        angleMoveOuter.setFixedCOM(false);
        integrator.getMoveManager().addMCMove(angleMoveOuter, 0.1);

        MCMoveClusterAngleMulti angleMoveMultiInner=null;
        if (armLength>2) {
            angleMoveMultiInner = new MCMoveClusterAngleMulti(pc, space, bonding, sim.getRandom(), 1,
                    new MCMoveClusterAngle.AtomChooserStarfl(sim.getRandom(), numArms, armLength, 2, (armLength+1) / 2), 5);
            angleMoveMultiInner.setBox(sim.box());
            angleMoveMultiInner.setFixedCOM(false);
            integrator.getMoveManager().addMCMove(angleMoveMultiInner, 0.4);
        }

        MCMoveClusterAngleMulti angleMoveMultiOuter = new MCMoveClusterAngleMulti(pc, space, bonding, sim.getRandom(), 1,
                new MCMoveClusterAngle.AtomChooserStarfl(sim.getRandom(), numArms, armLength, (armLength+1) / 2 + 1, armLength), 10);
        angleMoveMultiOuter.setBox(sim.box());
        angleMoveMultiOuter.setFixedCOM(false);
        integrator.getMoveManager().addMCMove(angleMoveMultiOuter, 0.4);

        // graphics block unchanged and still disabled
        if (false) {
            double size = (2*armLength + 5) * 1.5;
            sim.box().getBoundary().setBoxSize(Vector.of(new double[]{size, size, size}));
            sim.getController().addActivity(new ActivityIntegrate(integrator), Long.MAX_VALUE, 10);
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, "Star Single", 1);
            DisplayBox displayBox0 = simGraphic.getDisplayBox(sim.box());
            displayBox0.setPixelUnit(new Pixel(300.0 / size));
            displayBox0.setShowBoundary(false);
            ((DisplayBoxCanvasG3DSys) displayBox0.canvas).setBackgroundColor(Color.WHITE);
            DiameterHashByType diameters = new DiameterHashByType();
            diameters.setDiameter(typeCore, params.coreSigma);
            displayBox0.setDiameterHash(diameters);
            simGraphic.makeAndDisplayFrame();
            return;
        }

        long t1 = System.nanoTime();

        // equilibration — identical to VirialStarSingle
        ActivityIntegrate ai = new ActivityIntegrate(integrator, steps/10);
        sim.getController().runActivityBlocking(ai);

        long t2 = System.nanoTime();
        System.out.println("equilibration finished: "+(t2-t1)/1e9);
        if(angleMoveInner!=null) System.out.println("Angle move step size Inner   " + angleMoveInner.getStepSize());
        System.out.println("Angle move step size Outer  " + angleMoveOuter.getStepSize());
        if(angleMoveMultiInner!=null) System.out.println("Angle move step size MultiInner   " + angleMoveMultiInner.getStepSize());
        System.out.println("Angle move step size MultiOuter   " + angleMoveMultiOuter.getStepSize());

        integrator.getMoveManager().setEquilibrating(false);

        // ---------------------------------------------------------------
        // CHANGE 4: replaced MeterRadiusGyration + accRg + results.csv
        //           with MeterArmProfile + accProfile + profile.csv
        //
        // REMOVED (these 4 lines from VirialStarSingle):
        //   MeterRadiusGyration meterRg = new MeterRadiusGyration(sim.box());
        //   AccumulatorAverageFixed accRg = new AccumulatorAverageFixed(steps/1000);
        //   DataPumpListener pumpRg = new DataPumpListener(meterRg, accRg, 10);
        //   integrator.getEventManager().addListener(pumpRg);
        //
        // ADDED (same pattern, new meter):
        MeterArmProfile meterProfile = new MeterArmProfile(sim.box(), numArms, armLength);
        AccumulatorAverageFixed accProfile = new AccumulatorAverageFixed(steps / 1000);
        DataPumpListener pumpProfile = new DataPumpListener(meterProfile, accProfile, 10);
        integrator.getEventManager().addListener(pumpProfile);
        // ---------------------------------------------------------------

        // CHANGE 5: removed XYZWriter block from VirialStarSingle
        // (not needed for structural profile output)

        ai = new ActivityIntegrate(integrator, steps);
        sim.getController().runActivityBlocking(ai);

        System.out.println();
        if(angleMoveInner!=null) System.out.println("Angle move acceptance Inner   " + angleMoveInner.getTracker().acceptanceProbability());
        System.out.println("Angle move acceptance Outer  " + angleMoveOuter.getTracker().acceptanceProbability());
        if(angleMoveMultiInner!=null) System.out.println("Angle move acceptance MultiInner   " + angleMoveMultiInner.getTracker().acceptanceProbability());
        System.out.println("Angle move acceptance MultiOuter   " + angleMoveMultiOuter.getTracker().acceptanceProbability());
        System.out.println();

        // ---------------------------------------------------------------
        // CHANGE 6: replaced Rg print + results.csv write
        //           with profile print + profile.csv write
        //
        // REMOVED from VirialStarSingle:
        //   double avgRg = accRg.getData(accRg.AVERAGE).getValue(0);
        //   double errRg = accRg.getData(accRg.ERROR).getValue(0);
        //   double corRg = accRg.getData(accRg.BLOCK_CORRELATION).getValue(0);
        //   System.out.println("Rg2: "+avgRg ...);
        //   System.out.println("Rg: " + Math.sqrt(avgRg) ...);
        //   ... results.csv write block ...
        //
        // ADDED:
        IData avgProfile = accProfile.getData(accProfile.AVERAGE);
        IData errProfile = accProfile.getData(accProfile.ERROR);
        IData corProfile = accProfile.getData(accProfile.BLOCK_CORRELATION);

        System.out.println("Arm radial profile r(i): monomer index i (1=closest to core, armLength=tip)");
        System.out.println("i, r_avg, r_err, r_cor");
        for (int k = 1; k <= armLength; k++) {
            System.out.printf("%d, %.6f, %.6f, %.6f%n",
                    k,
                    avgProfile.getValue(k - 1),
                    errProfile.getValue(k - 1),
                    corProfile.getValue(k - 1));
        }

        // write profile.csv — one row per bead index per run
        // columns: numArms, armLength, coreSigma, epsCore, temperature, i, r_avg, r_err, r_cor
        String profileFile = System.getProperty("user.dir") + "/profile.csv";
        boolean fileExists = new java.io.File(profileFile).exists();
        try (FileWriter fw = new FileWriter(profileFile, true)) {
            if (!fileExists) {
                fw.write("numArms,armLength,coreSigma,epsCore,temperature,i,r_avg,r_err,r_cor\n");
            }
            for (int k = 1; k <= armLength; k++) {
                fw.write(String.format("%d,%d,%.2f,%.2f,%.4f,%d,%.8f,%.8f,%.8f%n",
                        numArms,
                        armLength,
                        params.coreSigma,
                        params.epsCore,
                        temperature,
                        k,
                        avgProfile.getValue(k - 1),
                        errProfile.getValue(k - 1),
                        corProfile.getValue(k - 1)));
            }
            System.out.println("Profile written to: " + profileFile);
        } catch (Exception e) {
            System.err.println("Warning: Could not write to profile.csv: " + e.getMessage());
        }
        // ---------------------------------------------------------------

        long t3 = System.nanoTime();
        System.out.println("time: "+(t3-t2)/1e9);
    }

    // VirialStarParams unchanged — identical to VirialStarSingle
    public static class VirialStarParams extends ParameterBase {
        public int armLength = 3;
        public double temperature = 1;
        public long numSteps = 1000000;
        public int numArms = 2;
        public boolean ideal = false;
        public double coreSigma = 1.0;
        public double epsCore   = 1.0;
        public boolean useCoreInHSRef = true;
        public double hsExtra = 0.0;
    }
}
