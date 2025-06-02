/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations.theta;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.DiameterHashByType;
import etomica.config.ConformationLinear;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.DataPumpListener;
import etomica.data.IData;
import etomica.data.meter.MeterRadiusGyration;
import etomica.data.types.DataGroup;
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
import etomica.space.BoundaryRectangularNonperiodic;
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
public class VirialChainSingleSQW {

    public static void main(String[] args) {
        VirialChainParams params = new VirialChainParams();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        } else {
            params.numSteps = 1000000;
            params.nH = 24;
            params.temperature = 3;
            params.lattice = true;
        }
        double temperature = params.temperature;
        long steps = params.numSteps;
        double lambda = params.lambda;
        int nH = params.nH;
        int nT = params.nT;
        if (nH < 0 || nT < 0) {
            throw new RuntimeException("nH: "+nH+"    and nT: "+nT);
        }
        int nSpheres = nH + nT;
        boolean doLattice = params.lattice;

        double epsH = 0.25;
        double epsT = 0.75;
        double epsHT = 0.25;
        System.out.println("epsH: "+epsH+"   epsT: "+epsT+"    epsHT: "+epsHT);

        Space space = Space3D.getInstance();

        ConformationLinear conf = new ConformationLinear(Space3D.getInstance(), 1, new double[]{0,Math.PI/2});
        AtomType typeH = AtomType.simple("H");
        AtomType typeT = AtomType.simple("T");
        ISpecies species = new SpeciesBuilder(Space3D.getInstance())
                .addCount(typeH, nH)
                .addCount(typeT, nT)
                .withConformation(conf).build();
        SpeciesManager sm = new SpeciesManager.Builder().addSpecies(species).build();

        double sigma = doLattice ? 0.9 : 1;
        lambda = doLattice ? 2 : lambda;

        PotentialMoleculePair pTarget = new PotentialMoleculePair(space, sm);
        System.out.println(nSpheres+"-mer chain");
        System.out.println("T: "+temperature);
        System.out.println("nH: "+nH);
        System.out.println("nT: "+nT);
        if (!doLattice) System.out.println("SQW lambda: "+lambda);

        P2HardGeneric p2HH = P2SquareWell.makePotential(sigma, lambda, epsH);
        P2HardGeneric p2TT = P2SquareWell.makePotential(sigma, lambda, epsH);
        P2HardGeneric p2HT = P2SquareWell.makePotential(sigma, lambda, epsHT);
        if (nH>0) pTarget.setAtomPotential(typeH, typeH, p2HH);
        if (nT>0) pTarget.setAtomPotential(typeT, typeT, p2TT);
        if (nH*nT>0) pTarget.setAtomPotential(typeH, typeT, p2HT);

        PotentialMasterBonding.FullBondingInfo bondingInfo = new PotentialMasterBonding.FullBondingInfo(sm);

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

        Simulation sim = new Simulation(Space3D.getInstance(), sm);
        sim.makeBox(new BoundaryRectangularNonperiodic(Space3D.getInstance()));
        sim.box().setNMolecules(species, 1);

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

        MCMoveClusterAngle angleMove = new MCMoveClusterAngle(pc, space, bonding, sim.getRandom(), 1);
        angleMove.setBox(sim.box());
        angleMove.setDoLattice(doLattice);
        integrator.getMoveManager().addMCMove(angleMove);

        MCMoveClusterReptate reptateMove = new MCMoveClusterReptate(pc, space, sim.getRandom());
        reptateMove.setBox(sim.box());
        reptateMove.setDoLattice(doLattice);
        integrator.getMoveManager().addMCMove(reptateMove);

        MCMoveClusterShuffle shuffleMove = new MCMoveClusterShuffle(pc, space, sim.getRandom());
        shuffleMove.setBox(sim.box());
        shuffleMove.setDoLattice(doLattice);
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
            DiameterHashByType diameterHash = (DiameterHashByType) displayBox0.getDiameterHash();
            diameterHash.setDiameter(typeH, 1);
            diameterHash.setDiameter(typeT, 1);

            simGraphic.makeAndDisplayFrame();

            return;
        }

        long t1 = System.nanoTime();
        // if running interactively, don't use the file
        ActivityIntegrate ai = new ActivityIntegrate(integrator, steps/10);
        sim.getController().runActivityBlocking(ai);

        System.out.println("equilibration finished");
        System.out.println("Angle move step size    " + angleMove.getStepSize());
        System.out.println("Shuffle move step size    "+shuffleMove.getStepSize());

        integrator.getMoveManager().setEquilibrating(false);

        DataSourceThetaSingle dsTheta = new DataSourceThetaSingle(sim.box(), temperature, 0, pcPair, null);
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
        System.out.println("Angle move acceptance " + angleMove.getTracker().acceptanceProbability());
        System.out.println("Shuffle move acceptance "+shuffleMove.getTracker().acceptanceProbability());
        System.out.println();

        DataGroup allYourBase = (DataGroup)accTheta.getData();
        IData avg = allYourBase.getData(accTheta.AVERAGE.index);
        IData err = allYourBase.getData(accTheta.ERROR.index);
        IData cor = allYourBase.getData(accTheta.BLOCK_CORRELATION.index);

        System.out.println("dlnqdb: "+avg.getValue(0)+"   err: "+err.getValue(0)+"   cor: "+cor.getValue(0));
        System.out.println("qbboq: "+avg.getValue(3)+"   err: "+err.getValue(3)+"   cor: "+cor.getValue(3));

        double avgRg = accRg.getData(accRg.AVERAGE).getValue(0);
        double errRg = accRg.getData(accRg.ERROR).getValue(0);
        double corRg = accRg.getData(accRg.BLOCK_CORRELATION).getValue(0);

        System.out.println("Rg2: "+avgRg+"   err: "+errRg+"  cor: "+corRg);

        long t2 = System.nanoTime();
        System.out.println("time: "+(t2-t1)/1e9);
    }

    /**
     * Inner class for parameters
     */
    public static class VirialChainParams extends ParameterBase {
        public double temperature = 1;
        public long numSteps = 1000000;
        public double lambda = 1.5;
        public int nH = 5;
        public int nT = 0;
        public boolean lattice;
    }
}
