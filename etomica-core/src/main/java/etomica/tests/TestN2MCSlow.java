/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.tests;

import etomica.action.BoxImposePbc;
import etomica.action.BoxInflate;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.DiameterHashByType;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.chem.elements.Nitrogen;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.DataPumpListener;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.meter.MeterPressure;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorBox;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveMolecule;
import etomica.integrator.mcmove.MCMoveRotateMolecule3D;
import etomica.lattice.LatticeCubicFcc;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.potential.*;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesBuilder;
import etomica.species.SpeciesGeneral;
import etomica.units.Electron;
import etomica.units.Kelvin;
import etomica.units.Pixel;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.util.random.RandomMersenneTwister;

import java.awt.*;

/**
 * Simple Lennard-Jones molecular dynamics simulation in 3D
 */
public class TestN2MCSlow extends Simulation {

    public PotentialMaster potentialMaster;
    public IntegratorMC integrator;
    public SpeciesGeneral species;
    public Box box;
    public MCMoveMolecule translateMove;
    public MCMoveRotateMolecule3D rotateMove;

    public TestN2MCSlow(int numMolecules, double temperatureK, boolean cellListing, Truncation trunc, boolean droplet) {
        super(Space3D.getInstance());
        setRandom(new RandomMersenneTwister(1));

        AtomType typeN = AtomType.element(Nitrogen.INSTANCE);
        AtomType typeM = AtomType.simple("M", 0);
        species = new SpeciesBuilder(space)
                .addAtom(typeN, Vector.of(-0.55, 0, 0))
                .addAtom(typeM, Vector.of(0, 0, 0))
                .addAtom(typeN, Vector.of(+0.55, 0, 0))
                .build();
        addSpecies(species);

        box = this.makeBox();

        double rc = 49.9999;
        potentialMaster = cellListing ? new PotentialMasterCell(this, rc) : new PotentialMaster();

        integrator = new IntegratorMC(this.getRandom(), potentialMaster, box);
        integrator.setTemperature(Kelvin.UNIT.toSim(temperatureK));

        translateMove = new MCMoveMolecule(potentialMaster, random, space, 0.2, 50);
        integrator.getMoveManager().addMCMove(translateMove);
        rotateMove = new MCMoveRotateMolecule3D(potentialMaster, random, space);
        integrator.getMoveManager().addMCMove(rotateMove);

        box.setNMolecules(species, numMolecules);
        new BoxInflate(box, space, 100.0 / (100 * 100 * 100)).actionPerformed();
        System.out.println("box size: " + box.getBoundary().getBoxSize());

        Potential2Soft pNNt, pNMt, pMMt;
        double sigma = 3.31;
        double epsilon = Kelvin.UNIT.toSim(36);
        double qN = Electron.UNIT.toSim(-0.482);
        double qM = -2 * qN;

        if (trunc == Truncation.EWALD) {
            throw new RuntimeException("You crazy");
        } else if (trunc == Truncation.CUT) {

            P2LennardJones pNNLJ = new P2LennardJones(space, sigma, epsilon);
            P2SoftSphericalTruncated pNNLJt = new P2SoftSphericalTruncated(space, pNNLJ, rc);
            P2SoftSphericalTruncatedForceShifted pNN1sf = new P2SoftSphericalTruncatedForceShifted(space, new P2Electrostatic(space, qN, qN), rc);

            P2Electrostatic pNM1 = new P2Electrostatic(space, qN, qM);
            pNMt = new P2SoftSphericalTruncatedForceShifted(space, pNM1, rc);

            P2Electrostatic pMM1 = new P2Electrostatic(space, qM, qM);
            pMMt = new P2SoftSphericalTruncatedForceShifted(space, pMM1, rc);

            potentialMaster.addPotential(pNNLJt, new AtomType[]{typeN, typeN});
            potentialMaster.addPotential(pNN1sf, new AtomType[]{typeN, typeN});
            potentialMaster.addPotential(pNMt, new AtomType[]{typeN, typeM});
            potentialMaster.addPotential(pMMt, new AtomType[]{typeM, typeM});
        } else {
            throw new RuntimeException("unknown truncation scheme");
        }


        if (!cellListing) {
            BoxImposePbc imposepbc = new BoxImposePbc(space);
            imposepbc.setBox(box);
            integrator.getEventManager().addListener(new IntegratorListenerAction(imposepbc, 10));
        }

        if (droplet) {
            double smallV = Math.PI * 4 / 3 * 11.3 * 11.3 * 11.3 * numMolecules / 100;
            Box initBox = makeBox(new BoundaryRectangularPeriodic(space, Math.pow(smallV, 1.0 / 3.0)));
            initBox.setNMolecules(species, numMolecules);

            ConfigurationLattice configuration = new ConfigurationLattice(new LatticeCubicFcc(space), space);
            configuration.initializeCoordinates(initBox);
            for (int i = 0; i < numMolecules; i++) {
                Vector r0 = initBox.getLeafList().get(3 * i).getPosition();
                for (int j = 0; j < 3; j++) {
                    IAtom initAtom = initBox.getLeafList().get(3 * i + j);
                    Vector dr = Vector.d(3);
                    dr.Ev1Mv2(initAtom.getPosition(), r0);
                    initBox.getBoundary().nearestImage(dr);
                    IAtom realAtom = box.getLeafList().get(3 * i + j);
                    realAtom.getPosition().Ev1Pv2(r0, dr);
                }
            }
            removeBox(initBox);
        } else {
            ConfigurationLattice configuration = new ConfigurationLattice(new LatticeCubicFcc(space), space);
            configuration.initializeCoordinates(box);
        }
    }

    public IntegratorBox getIntegrator() {
        return integrator;
    }

    public static void main(String[] args) {
        TestN2MCSlow.SimParams params = new TestN2MCSlow.SimParams();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        } else {
            // stuff
        }
        int numMolecules = params.numMolecules;

        final TestN2MCSlow sim = new TestN2MCSlow(numMolecules, params.temperatureK, false, params.trunc, params.droplet);

        boolean graphic = false;
        if (graphic) {
            sim.getController().addActivity(new ActivityIntegrate(sim.integrator));
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, "N2");
            DiameterHashByType diameters = (DiameterHashByType) simGraphic.getDisplayBox(sim.box).getDiameterHash();
            AtomType typeN = sim.species.getTypeByName("N");
            diameters.setDiameter(typeN, 3.31);
            ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box).getColorScheme()).setColor(typeN, Color.BLUE);
            simGraphic.getDisplayBox(sim.box).setPixelUnit(new Pixel(4));
            simGraphic.makeAndDisplayFrame();
            return;
        }

//        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, params.numSteps / 10));
//        sim.integrator.resetStepCount();

        int bs = params.numSteps / (100 * 2 * numMolecules);
        if (bs == 0) bs = 1;
        MeterPressure pMeter = new MeterPressure(sim.space);
        pMeter.setIntegrator(sim.integrator);
        AccumulatorAverage pAccumulator = new AccumulatorAverageFixed(bs);
        DataPumpListener pPump = new DataPumpListener(pMeter, pAccumulator, 2 * numMolecules);
        sim.integrator.getEventManager().addListener(pPump);

        bs = params.numSteps / (100 * 10);
        if (bs == 0) bs = 1;
        MeterPotentialEnergyFromIntegrator energyMeter = new MeterPotentialEnergyFromIntegrator(sim.integrator);
        AccumulatorAverage energyAccumulator = new AccumulatorAverageFixed(bs);
        DataPumpListener energyPump = new DataPumpListener(energyMeter, energyAccumulator, 10);
        sim.integrator.getEventManager().addListener(energyPump);

        long t1 = System.nanoTime();
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, params.numSteps));
        long t2 = System.nanoTime();
        System.out.println("time: " + (t2 - t1) / 1e9);

        System.out.println("Translate step size: " + sim.translateMove.getStepSize());
        System.out.println("Rotate step size: " + sim.rotateMove.getStepSize());

        double avgP = pAccumulator.getData(pAccumulator.AVERAGE).getValue(0);
        double errP = pAccumulator.getData(pAccumulator.ERROR).getValue(0);
        double corP = pAccumulator.getData(pAccumulator.BLOCK_CORRELATION).getValue(0);
        System.out.println("P " + avgP + " " + errP + " " + corP);

        double avgPE = energyAccumulator.getData(energyAccumulator.AVERAGE).getValue(0) / numMolecules;
        double errPE = energyAccumulator.getData(energyAccumulator.ERROR).getValue(0) / numMolecules;
        double corPE = energyAccumulator.getData(energyAccumulator.BLOCK_CORRELATION).getValue(0);
        System.out.println("PE " + avgPE + " " + errPE + " " + corPE);
    }

    public enum Truncation {
        CUT, EWALD
    }

    public static class SimParams extends ParameterBase {
        public int numMolecules = 100;
        public int numSteps = 400000;
        public double temperatureK = 50;
        public boolean droplet = true;
        public Truncation trunc = Truncation.CUT;
    }
}
