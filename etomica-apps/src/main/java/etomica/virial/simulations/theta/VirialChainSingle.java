/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations.theta;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.config.ConformationLinear;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.DataPumpListener;
import etomica.data.IData;
import etomica.data.meter.MeterRadiusGyration;
import etomica.data.types.DataGroup;
import etomica.graphics.ColorSchemeRandom;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayBoxCanvasG3DSys;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveStepTracker;
import etomica.potential.*;
import etomica.potential.compute.NeighborManagerIntra;
import etomica.potential.compute.PotentialComputeAggregate;
import etomica.potential.compute.PotentialComputePair;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesBuilder;
import etomica.species.SpeciesManager;
import etomica.units.Pixel;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.util.collections.IntArrayList;
import etomica.virial.mcmove.MCMoveClusterAngle;

import java.awt.*;
import java.util.ArrayList;
import java.util.List;

/**
 * Single-molecule calculations for input into theta temperature calculations
 */
public class VirialChainSingle {

    public static void main(String[] args) {
        VirialChainParams params = new VirialChainParams();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        } else {
            params.nSpheres = 128;
            params.temperature = 4; //.42304894079030;
            params.kBend = 0;
            params.numSteps = 1000000;
        }
        int nSpheres = params.nSpheres;
        double temperature = params.temperature;
        long steps = params.numSteps;
        double eFENE = params.eFENE;
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
        System.out.println(nSpheres+"-mer chain at T = "+temperature);
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
        pTarget.setAtomPotential(type, type, p2);

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
        if (kBend < Double.POSITIVE_INFINITY) {
            P3BondAngleStiffChain p3 = new P3BondAngleStiffChain(kBend);
            bondingInfo.setBondingPotentialTriplet(species, p3, triplets);
        }

        Simulation sim = new Simulation(Space3D.getInstance(), sm);
        sim.makeBox();
        sim.box().setNMolecules(species, 1);

        PotentialMasterBonding.FullBondingInfo bondingInfodk = new PotentialMasterBonding.FullBondingInfo(sm);

        PotentialMasterBonding pmBondingdk = new PotentialMasterBonding(sm, sim.box(), bondingInfodk);
        P3BondAngleStiffChain p3dk = new P3BondAngleStiffChain(1);
        bondingInfodk.setBondingPotentialTriplet(species, p3dk, triplets);

        PotentialMasterBonding pmBonding = new PotentialMasterBonding(sm, sim.box(), bondingInfo);
        PotentialComputePair pcPair = new PotentialComputePair(sm, sim.box(), new NeighborManagerIntra(sim.box(), bondingInfo), pTarget.getAtomPotentials());
        PotentialComputeAggregate pc = new PotentialComputeAggregate(pmBonding, pcPair);

        IntegratorMC integrator = new IntegratorMC(pc, sim.getRandom(), temperature, sim.box());

        System.out.println(steps+" steps");

        IntArrayList[] bonding = new IntArrayList[nSpheres];
        bonding[0] = new IntArrayList(new int[]{1});
        for (int i=1; i<nSpheres-1; i++) {
            bonding[i] = new IntArrayList(new int[]{i-1,i+1});
        }
        bonding[nSpheres-1] = new IntArrayList(new int[]{nSpheres-2});

        MCMoveClusterAngle angleMove1 = null, angleMove2 = null, angleMove3 = null, angleMove4 = null;
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

        MCMoveClusterReptate reptateMove = new MCMoveClusterReptate(pc, space, sim.getRandom());
        reptateMove.setBox(sim.box());
        integrator.getMoveManager().addMCMove(reptateMove);

        MCMoveClusterShuffle shuffleMove = new MCMoveClusterShuffle(pc, space, sim.getRandom());
        shuffleMove.setBox(sim.box());
        integrator.getMoveManager().addMCMove(shuffleMove);
        ((MCMoveStepTracker)shuffleMove.getTracker()).setAcceptanceTarget(0.3);

        if (false) {
            ActivityIntegrate ai = new ActivityIntegrate(integrator);
            sim.getController().addActivity(ai, Long.MAX_VALUE, 10);

            double size = (nSpheres + 5) * 1.5;
            sim.box().getBoundary().setBoxSize(Vector.of(new double[]{size, size, size}));
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, "stuff", 1);
            DisplayBox displayBox0 = simGraphic.getDisplayBox(sim.box());
            displayBox0.setPixelUnit(new Pixel(300.0 / size));
            displayBox0.setShowBoundary(false);
            ((DisplayBoxCanvasG3DSys) displayBox0.canvas).setBackgroundColor(Color.WHITE);

            ColorSchemeRandom colorScheme = new ColorSchemeRandom(sim.box(), sim.getRandom());
            displayBox0.setColorScheme(colorScheme);

            simGraphic.makeAndDisplayFrame();

            return;
        }

        long t1 = System.nanoTime();
        // if running interactively, don't use the file
        ActivityIntegrate ai = new ActivityIntegrate(integrator, steps/10);
        sim.getController().runActivityBlocking(ai);

        System.out.println("equilibration finished");
        if (angleMove2 == null) {
            System.out.println("Angle move step size    " + angleMove1.getStepSize());
        }
        else if (angleMove3 == null) {
            System.out.println("Angle move step size    " + angleMove1.getStepSize() + " " + angleMove2.getStepSize());
        }
        else {
            System.out.println("Angle move step size    " + angleMove1.getStepSize() + " " + angleMove2.getStepSize() + " " + angleMove3.getStepSize() + " " + angleMove4.getStepSize());
        }
        System.out.println("Shuffle move step size    "+shuffleMove.getStepSize());

        integrator.getMoveManager().setEquilibrating(false);

        DataSourceThetaSingle dsTheta = new DataSourceThetaSingle(sim.box(), temperature, kBend, pcPair, pmBondingdk);
        AccumulatorAverageFixed accTheta = new AccumulatorAverageFixed(steps/1000);
        DataPumpListener pumpTheta = new DataPumpListener(dsTheta, accTheta, 10);
        integrator.getEventManager().addListener(pumpTheta);

        MeterRadiusGyration meterRg = new MeterRadiusGyration(sim.box());
        AccumulatorAverageFixed accRg = new AccumulatorAverageFixed(steps/1000);
        DataPumpListener pumpRg = new DataPumpListener(meterRg, accRg, 10);
        integrator.getEventManager().addListener(pumpRg);

        ai = new ActivityIntegrate(integrator, steps);
        sim.getController().runActivityBlocking(ai);

        System.out.println();
        System.out.println("Reptate move acceptance "+reptateMove.getTracker().acceptanceProbability());
        if (angleMove2 == null) {
            System.out.println("Angle move acceptance " + angleMove1.getTracker().acceptanceProbability());
        }
        else if (angleMove3 == null) {
            System.out.println("Angle move acceptance " + angleMove1.getTracker().acceptanceProbability() + " " + angleMove2.getTracker().acceptanceProbability());
        }
        else {
            System.out.println("Angle move acceptance " + angleMove1.getTracker().acceptanceProbability() + " " + angleMove2.getTracker().acceptanceProbability() + " " + angleMove3.getTracker().acceptanceProbability() + " " + angleMove4.getTracker().acceptanceProbability());
        }
        System.out.println("Shuffle move acceptance "+shuffleMove.getTracker().acceptanceProbability());
        System.out.println();

        DataGroup allYourBase = (DataGroup)accTheta.getData();
        IData avg = allYourBase.getData(accTheta.AVERAGE.index);
        IData err = allYourBase.getData(accTheta.ERROR.index);
        IData cor = allYourBase.getData(accTheta.BLOCK_CORRELATION.index);

        System.out.println("dlnqdb: "+avg.getValue(0)+"   err: "+err.getValue(0)+"   cor: "+cor.getValue(0));
        System.out.println("dlnqdk: "+avg.getValue(1)+"   err: "+err.getValue(1)+"   cor: "+cor.getValue(1));
        System.out.println("qkboq: "+avg.getValue(2)+"   err: "+err.getValue(2)+"   cor: "+cor.getValue(2));
        System.out.println("qbboq: "+avg.getValue(3)+"   err: "+err.getValue(3)+"   cor: "+cor.getValue(3));

        double avgRg = accRg.getData(accRg.AVERAGE).getValue(0);
        double errRg = accRg.getData(accRg.ERROR).getValue(0);
        double corRg = accRg.getData(accRg.BLOCK_CORRELATION).getValue(0);

        System.out.println("Rg2: "+avgRg+"   err: "+errRg+"  cor: "+corRg);

        long t2 = System.nanoTime();
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
        public double eFENE = 0;
        public double rc = 0;
        public double bondLength = 1;
        public TruncationChoice truncation = TruncationChoice.SHIFT;
        public double kBend = 0;
    }
}
