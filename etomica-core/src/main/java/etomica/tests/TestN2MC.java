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
import etomica.data.meter.MeterPotentialEnergyFromIntegratorFasterer;
import etomica.data.meter.MeterPressureFasterer;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorBoxFasterer;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorMCFasterer;
import etomica.integrator.mcmove.MCMoveMoleculeFasterer;
import etomica.integrator.mcmove.MCMoveMoleculeRotateFasterer;
import etomica.lattice.LatticeCubicFcc;
import etomica.nbr.cell.PotentialMasterCellFasterer;
import etomica.potential.*;
import etomica.potential.compute.PotentialComputeAggregate;
import etomica.potential.compute.PotentialComputeEwaldFourier;
import etomica.potential.ewald.P2Ewald1Real;
import etomica.potential.ewald.P2Ewald6Real;
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
 * Monte Carlo simulation for TraPPE N2.  Can use Ewald or truncated potential;
 * can initialize on a lattice filling the box or one at a liquid density in
 * middle (forms a droplet).  Parameters can be set to match the n2_droplet
 * example from Cassandra.
 */
public class TestN2MC extends Simulation {

    public PotentialComputeEwaldFourier ewald;
    public PotentialMasterFasterer potentialMaster;
    public PotentialComputeAggregate pcAggregate;
    public IntegratorMCFasterer integrator;
    public SpeciesGeneral species;
    public Box box;
    public MCMoveMoleculeFasterer translateMove;
    public MCMoveMoleculeRotateFasterer rotateMove;

    public TestN2MC(int numMolecules, double temperatureK, boolean cellListing, Truncation trunc, boolean droplet) {
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

        potentialMaster = cellListing ? new PotentialMasterCellFasterer(this, box, 3, BondingInfo.noBonding()) : new PotentialMasterFasterer(this, box, BondingInfo.noBonding());

        boolean doEwald = (trunc == Truncation.EWALD || trunc == Truncation.EWALD1);
        ewald = doEwald ? new PotentialComputeEwaldFourier(this, box, BondingInfo.noBonding()) : null;

        if (doEwald) {
            pcAggregate = new PotentialComputeAggregate(potentialMaster, ewald);
        } else {
            pcAggregate = new PotentialComputeAggregate(potentialMaster);
        }

        integrator = new IntegratorMCFasterer(this, pcAggregate, box);
        integrator.setTemperature(Kelvin.UNIT.toSim(temperatureK));

        translateMove = new MCMoveMoleculeFasterer(random, pcAggregate, box);
        translateMove.setStepSizeMax(50);
        integrator.getMoveManager().addMCMove(translateMove);
        rotateMove = new MCMoveMoleculeRotateFasterer(random, pcAggregate, box);
        integrator.getMoveManager().addMCMove(rotateMove);

        box.setNMolecules(species, numMolecules);
        new BoxInflate(box, space, 100.0 / (100 * 100 * 100)).actionPerformed();
        System.out.println("box size: " + box.getBoundary().getBoxSize());

        Potential2Soft pNNt, pNMt, pMMt;
        double sigma = 3.31;
        double epsilon = Kelvin.UNIT.toSim(36);
        double qN = Electron.UNIT.toSim(-0.482);
        double qM = -2 * qN;

        if (doEwald) {
            // 3.30703 gives rc=50
            // Cassandra with precision=1e-5 gives s~=3.4
            double s = 3;
            PotentialComputeEwaldFourier.EwaldParams p = ewald.getOptimalParams(s, 4);
            System.out.println("optimal: " + p);

//            p.rCut = 49.9999;
//            p.alpha = s / p.rCut;
//            p.kCut = 2 * p.alpha * s;
            System.out.print("rc: " + p.rCut + "\n" +
                    "kc: " + p.kCut + "\n" +
                    "alpha: " + p.alpha + "\n");

            ewald.setkCut(p.kCut);
            ewald.setAlpha(p.alpha);
            ewald.setCharge(typeN, qN);
            ewald.setCharge(typeM, qM);

            P2Ewald1Real pNN1 = new P2Ewald1Real(qN * qN, p.alpha);
            if (trunc == Truncation.EWALD) {
                ewald.setAlpha6(p.alpha);
                ewald.setR6Coefficient(typeN, sigma, epsilon);

                P2SoftSphere pNN12 = new P2SoftSphere(space, sigma, 4 * epsilon, 12);
                P2Ewald6Real pNN6 = new P2Ewald6Real(sigma, epsilon, sigma, epsilon, p.alpha);
                pNNt = new P2SoftSphericalSumTruncated(space, p.rCut, new Potential2SoftSpherical[]{pNN12, pNN6, pNN1});
            } else {
                P2LennardJones pNNLJ = new P2LennardJones(space, sigma, epsilon);
                pNNt = new P2SoftSphericalSumTruncated(space, p.rCut, pNNLJ, pNN1);
            }

            P2Ewald1Real pNM1 = new P2Ewald1Real(qN * qM, p.alpha);
            pNMt = new P2SoftSphericalTruncated(space, pNM1, p.rCut);

            P2Ewald1Real pMM1 = new P2Ewald1Real(qM * qM, p.alpha);
            pMMt = new P2SoftSphericalTruncated(space, pMM1, p.rCut);
        } else if (trunc == Truncation.CUT) {

            double rc = 49.9999;

            if (!droplet) {
                // let's pretend force-shifting is adequate
                P2LennardJones pNNLJ = new P2LennardJones(space, sigma, epsilon);
                P2SoftSphericalTruncatedForceShifted pNN1sf = new P2SoftSphericalTruncatedForceShifted(space, new P2Electrostatic(space, qN, qN), rc);

                pNNt = new P2SoftSphericalSumTruncated(space, rc, pNNLJ, pNN1sf);

                P2Electrostatic pNM1 = new P2Electrostatic(space, qN, qM);
                pNMt = new P2SoftSphericalSumTruncatedForceShifted(space, rc, pNM1);

                P2Electrostatic pMM1 = new P2Electrostatic(space, qM, qM);
                pMMt = new P2SoftSphericalSumTruncatedForceShifted(space, rc, pMM1);
            } else {
                // don't even need to worry about force-shifting because molecules can't escape the droplet
                P2LennardJones pNNLJ = new P2LennardJones(space, sigma, epsilon);
                P2SoftSphericalSumTruncated pNN1sf = new P2SoftSphericalSumTruncated(space, rc, new P2Electrostatic(space, qN, qN));

                pNNt = new P2SoftSphericalSumTruncated(space, rc, pNNLJ, pNN1sf);

                P2Electrostatic pNM1 = new P2Electrostatic(space, qN, qM);
                pNMt = new P2SoftSphericalSumTruncated(space, rc, pNM1);

                P2Electrostatic pMM1 = new P2Electrostatic(space, qM, qM);
                pMMt = new P2SoftSphericalSumTruncated(space, rc, pMM1);
            }
        } else {
            throw new RuntimeException("unknown truncation scheme");
        }

        potentialMaster.setPairPotential(typeN, typeN, pNNt);
        potentialMaster.setPairPotential(typeN, typeM, pNMt);
        potentialMaster.setPairPotential(typeM, typeM, pMMt);

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

    public IntegratorBoxFasterer getIntegrator() {
        return integrator;
    }

    public static void main(String[] args) {
        TestN2MC.SimParams params = new TestN2MC.SimParams();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        } else {
//            params.trunc = Truncation.EWALD1;
        }
        int numMolecules = params.numMolecules;

        final TestN2MC sim = new TestN2MC(numMolecules, params.temperatureK, false, params.trunc, params.droplet);

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

        if (false && !params.droplet) {
            sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, params.numSteps / 10));
            sim.integrator.resetStepCount();
            System.out.println("equilibration done");
        }

        int bs = params.numSteps / (100 * 2 * numMolecules);
        if (bs == 0) bs = 1;
        MeterPressureFasterer pMeter = new MeterPressureFasterer(sim.box, sim.pcAggregate);
        pMeter.setTemperature(sim.integrator.getTemperature());
        AccumulatorAverage pAccumulator = new AccumulatorAverageFixed(bs);
        DataPumpListener pPump = new DataPumpListener(pMeter, pAccumulator, 2 * numMolecules);
        sim.integrator.getEventManager().addListener(pPump);

        bs = params.numSteps / (100 * 10);
        if (bs == 0) bs = 1;
//        MeterPotentialEnergyFasterer energyMeter = new MeterPotentialEnergyFasterer(sim.box, sim.pcAggregate);
        MeterPotentialEnergyFromIntegratorFasterer energyMeter = new MeterPotentialEnergyFromIntegratorFasterer(sim.integrator);
        AccumulatorAverage energyAccumulator = new AccumulatorAverageFixed(bs);
        DataPumpListener energyPump = new DataPumpListener(energyMeter, energyAccumulator, 10);
        sim.integrator.getEventManager().addListener(energyPump);

        long t1 = System.nanoTime();
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, params.numSteps));
        long t2 = System.nanoTime();
        System.out.println("time: " + (t2 - t1) / 1e9);
        if (sim.ewald != null)
            System.out.println("Fourier: " + sim.ewald.fNum + " " + (sim.ewald.fTime / (double) sim.ewald.fNum));
        System.out.println("Pair: " + sim.potentialMaster.numAll + " " + (sim.potentialMaster.tAll / (double) sim.potentialMaster.numAll));

        if (sim.ewald != null)
            System.out.println("Fourier one: " + sim.ewald.numMC + " " + (sim.ewald.tMC / (double) sim.ewald.numMC));
        System.out.println("Pair one: " + sim.potentialMaster.numMC + " " + (sim.potentialMaster.tMC / (double) sim.potentialMaster.numMC));

        System.out.println("Translate step size: " + sim.translateMove.getStepSize() + " " + sim.translateMove.getTracker().acceptanceRatio());
        System.out.println("Rotate step size: " + sim.rotateMove.getStepSize() + " " + sim.rotateMove.getTracker().acceptanceRatio());

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
        CUT, EWALD, EWALD1
    }

    public static class SimParams extends ParameterBase {
        // these parameters correspond to Cassandra's n2_droplet example
        public int numMolecules = 100;
        public int numSteps = 400000;
        public double temperatureK = 50;
        public boolean droplet = true;
        public Truncation trunc = Truncation.CUT;
    }
}
