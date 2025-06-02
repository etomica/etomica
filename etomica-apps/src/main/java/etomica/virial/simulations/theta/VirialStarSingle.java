/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations.theta;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.DataPumpListener;
import etomica.data.DataSourceCountSteps;
import etomica.data.meter.MeterRadiusGyration;
import etomica.graphics.*;
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
import etomica.virial.mcmove.MCMoveClusterStretch;

import java.awt.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;


public class VirialStarSingle {

    public static void main(String[] args) {
        VirialStarParams params = new VirialStarParams();

        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        } else {
            params.armLength = 4;
            params.numArms = 4;
            params.temperature = 5;
            params.numSteps = 1000000;
        }
        int numArms = params.numArms;
        int armLength = params.armLength;
        double temperature = params.temperature;
        long steps = params.numSteps;
        boolean ideal = params.ideal;

        Space space = Space3D.getInstance();

        ConformationStarPolymerGraft conf = new ConformationStarPolymerGraft(Space3D.getInstance(), numArms, armLength);
        conf.setSigma0(1);
        AtomType type = AtomType.simple("A");
        ISpecies species = new SpeciesBuilder(Space3D.getInstance())
                .addCount(type, 1+numArms*armLength)
                .withConformation(conf).build();
        SpeciesManager sm = new SpeciesManager.Builder().addSpecies(species).build();

        PotentialMoleculePair pTarget = new PotentialMoleculePair(space, sm);
        System.out.println(numArms+" arms of length "+armLength+" at T = "+temperature);
        IPotential2 p2 = new P2LennardJones(1, 1);

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


        Simulation sim = new Simulation(Space3D.getInstance(), sm);
        sim.makeBox(new BoundaryRectangularNonperiodic(space));
        sim.box().setNMolecules(species, 1);

        System.out.println("random seeds: "+ Arrays.toString(sim.getRandomSeeds()));

        pTarget.setAtomPotential(type, type, p2);

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

        MCMoveClusterAngle angleMove = new MCMoveClusterAngle(pc, space, bonding, sim.getRandom(), 1);
        angleMove.setBox(sim.box());
        integrator.getMoveManager().addMCMove(angleMove);

        MCMoveClusterShuffle shuffleMove = null;
        if (armLength > 5) {
            shuffleMove = new MCMoveClusterShuffle(pc, space, sim.getRandom());
            shuffleMove.setBox(sim.box());
            shuffleMove.setBonding(bonding);
            shuffleMove.setStepSizeMax(armLength - 2);
            integrator.getMoveManager().addMCMove(shuffleMove);
            ((MCMoveStepTracker) shuffleMove.getTracker()).setAcceptanceTarget(0.3);
        }

        if (false) {
            double size = (2*armLength + 5) * 1.5;
            sim.box().getBoundary().setBoxSize(Vector.of(new double[]{size, size, size}));
            sim.getController().addActivity(new ActivityIntegrate(integrator), Long.MAX_VALUE, 10);
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, "Star Single", 1);
            DisplayBox displayBox0 = simGraphic.getDisplayBox(sim.box());
            displayBox0.setPixelUnit(new Pixel(300.0 / size));
            displayBox0.setShowBoundary(false);
            ((DisplayBoxCanvasG3DSys) displayBox0.canvas).setBackgroundColor(Color.WHITE);
            ColorScheme colorScheme = new ColorScheme() {
                @Override
                public Color getAtomColor(IAtom a) {
                    Color[] c = new Color[]{Color.BLACK, Color.RED, Color.BLUE, Color.GREEN, Color.YELLOW, Color.CYAN, Color.MAGENTA, Color.ORANGE, new Color(128,128,128)};
                    return c[(a.getIndex()+armLength-1)/armLength];
                }
            };
            displayBox0.setColorScheme(colorScheme);

            simGraphic.makeAndDisplayFrame();

            final DisplayTextBox stepsBox = new DisplayTextBox();
            stepsBox.setLabel("Steps");
            DataSourceCountSteps stepCounter = new DataSourceCountSteps(integrator);
            DataPumpListener stepPump = new DataPumpListener(stepCounter, stepsBox, 1000);
            integrator.getEventManager().addListener(stepPump);
            simGraphic.add(stepsBox);

            return;
        }


        long t1 = System.nanoTime();

        // if running interactively, don't use the file
        ActivityIntegrate ai = new ActivityIntegrate(integrator, steps/10);
        sim.getController().runActivityBlocking(ai);

        long t2 = System.nanoTime();
        System.out.println("equilibration finished: "+(t2-t1)/1e9);
        System.out.println("Angle move step size    " + angleMove.getStepSize());
        if (shuffleMove!=null) System.out.println("Shuffle move step size    "+shuffleMove.getStepSize());

        integrator.getMoveManager().setEquilibrating(false);

        MeterRadiusGyration meterRg = new MeterRadiusGyration(sim.box());
        AccumulatorAverageFixed accRg = new AccumulatorAverageFixed(steps/1000);
        DataPumpListener pumpRg = new DataPumpListener(meterRg, accRg, 10);
        integrator.getEventManager().addListener(pumpRg);

        ai = new ActivityIntegrate(integrator, steps);
        sim.getController().runActivityBlocking(ai);

        System.out.println();

        System.out.println("Angle move acceptance " + angleMove.getTracker().acceptanceProbability());
        if (shuffleMove!=null) System.out.println("Shuffle move acceptance " + shuffleMove.getTracker().acceptanceProbability());
        System.out.println();

        double avgRg = accRg.getData(accRg.AVERAGE).getValue(0);
        double errRg = accRg.getData(accRg.ERROR).getValue(0);
        double corRg = accRg.getData(accRg.BLOCK_CORRELATION).getValue(0);

        System.out.println("Rg2: "+avgRg+"   err: "+errRg+"  cor: "+corRg);

        long t3 = System.nanoTime();

        System.out.println("time: "+(t3-t2)/1e9);
    }

    /**
     * Inner class for parameters
     */
    public static class VirialStarParams extends ParameterBase {
        public int armLength = 3;
        public double temperature = 1;
        public long numSteps = 1000000;
        public int numArms = 2;
        public boolean ideal = false;
    }
}
