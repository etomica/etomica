/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.tests;

import etomica.action.BoxInflate;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.DiameterHashByType;
import etomica.box.Box;
import etomica.chem.elements.Copper;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.DataPumpListener;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.meter.MeterPressure;
import etomica.data.meter.MeterPressureFromIntegrator;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.lattice.LatticeCubicFcc;
import etomica.nbr.list.NeighborListManager;
import etomica.potential.BondingInfo;
import etomica.potential.EmbeddingSqrt;
import etomica.potential.P2LREPPhi;
import etomica.potential.P2LREPV;
import etomica.potential.compute.PotentialComputeEAM;
import etomica.simulation.Simulation;
import etomica.space3d.Space3D;
import etomica.species.SpeciesGeneral;
import etomica.units.ElectronVolt;
import etomica.units.Kelvin;
import etomica.units.Pascal;
import etomica.units.Pixel;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

/**
 * LREP Cu simulation with an FCC crystal.
 */
public class TestLREP extends Simulation {

    public IntegratorVelocityVerlet integrator;
    public SpeciesGeneral species;
    public Box box;

    public TestLREP(int numAtoms, double temperature, double density) {
        super(Space3D.getInstance());

        species = SpeciesGeneral.monatomic(space, AtomType.element(Copper.INSTANCE), true);
        addSpecies(species);

        box = this.makeBox();
        NeighborListManager nbrs = new NeighborListManager(this.getSpeciesManager(), box, 2, 7.2, BondingInfo.noBonding());
        nbrs.setDoDownNeighbors(true);
        PotentialComputeEAM potentialMaster = new PotentialComputeEAM(getSpeciesManager(), box, nbrs);
        potentialMaster.doAllTruncationCorrection = false;
        integrator = new IntegratorVelocityVerlet(potentialMaster, random, 0.001, temperature, box);
        integrator.setIsothermal(true);

        box.setNMolecules(species, numAtoms);
        BoxInflate inflater = new BoxInflate(box, space);
        inflater.setTargetDensity(density);
        inflater.actionPerformed();

        AtomType leafType = species.getLeafType();
        P2LREPV p2 = new P2LREPV();
        potentialMaster.setPairPotential(leafType, leafType, p2);
        P2LREPPhi pRho = new P2LREPPhi();
        potentialMaster.setRhoPotential(leafType, pRho);
        EmbeddingSqrt f = new EmbeddingSqrt(1);
        potentialMaster.setEmbeddingPotential(leafType, f);

        new ConfigurationLattice(new LatticeCubicFcc(space), space).initializeCoordinates(box);
    }

    public static void main(String[] args) {
        SimParams params = new SimParams();
        ParseArgs.doParseArgs(params, args);
        int numAtoms = params.numAtoms;
        double density = params.density;
        double temperature = Kelvin.UNIT.toSim(params.temperatureK);

        TestLREP sim = new TestLREP(numAtoms, temperature, density);
        sim.integrator.reset();

        double uLat = ElectronVolt.UNIT.fromSim(sim.integrator.getPotentialEnergy()/numAtoms);
        System.out.println("uLat (eV/atom): " + uLat);

        MeterPressure pMeterLat = new MeterPressure(sim.box, sim.integrator.getPotentialCompute());
        pMeterLat.setTemperature(0); // Lattice
        double pLat = 1.0E-9*Pascal.UNIT.fromSim(pMeterLat.getDataAsScalar());
        System.out.println("pLat (GPa): " + pLat);




        //??????????
        sim.integrator.reset();




        int steps = params.numSteps / numAtoms;
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, steps / 10));

        int bs = steps / (4 * 100);
        MeterPressureFromIntegrator pMeter = new MeterPressureFromIntegrator(sim.integrator);
        AccumulatorAverage pAccumulator = new AccumulatorAverageFixed(bs);
        DataPumpListener pPump = new DataPumpListener(pMeter, pAccumulator, 4);
        sim.integrator.getEventManager().addListener(pPump);

        MeterPotentialEnergyFromIntegrator energyMeter = new MeterPotentialEnergyFromIntegrator(sim.integrator);
        AccumulatorAverage energyAccumulator = new AccumulatorAverageFixed(bs);
        DataPumpListener energyPump = new DataPumpListener(energyMeter, energyAccumulator, 4);
        sim.integrator.getEventManager().addListener(energyPump);


        boolean graphic = false;
        if (graphic) {
            sim.getController().addActivity(new ActivityIntegrate(sim.integrator));
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, "Cu");
            DiameterHashByType diameters = (DiameterHashByType) simGraphic.getDisplayBox(sim.box).getDiameterHash();
            AtomType typeN = sim.species.getTypeByName("Cu");
            diameters.setDiameter(typeN, 1.5);
//            ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box).getColorScheme()).setColor(typeN, Color.BLUE);
            simGraphic.getDisplayBox(sim.box).setPixelUnit(new Pixel(4));
            simGraphic.makeAndDisplayFrame();
            return;
        }




        long t1 = System.currentTimeMillis();
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, steps));
        long t2 = System.currentTimeMillis();

        System.out.println("runtime: " + (t2 - t1) * 0.001);

        double avgP = 1E-9*Pascal.UNIT.fromSim(pAccumulator.getData(pAccumulator.AVERAGE).getValue(0));
        double errP = 1E-9*Pascal.UNIT.fromSim(pAccumulator.getData(pAccumulator.ERROR).getValue(0));
        double corP = pAccumulator.getData(pAccumulator.BLOCK_CORRELATION).getValue(0);
        System.out.println("P (GPa): " + avgP + " " + errP + " " + corP);

        double avgPE = ElectronVolt.UNIT.fromSim(energyAccumulator.getData(energyAccumulator.AVERAGE).getValue(0)/numAtoms);
        double errPE = ElectronVolt.UNIT.fromSim(energyAccumulator.getData(energyAccumulator.ERROR).getValue(0)/numAtoms);
        double corPE = energyAccumulator.getData(energyAccumulator.BLOCK_CORRELATION).getValue(0);
        System.out.println("PE (eV/atom): " + avgPE + " " + errPE + " " + corPE);

        double uHarm = ElectronVolt.UNIT.fromSim(1.5*(numAtoms-1.0)*temperature/numAtoms);
        System.out.println("uAnhC(eV/atom): " + (avgPE-uLat-uHarm) + " " + errPE + " " + corPE);

        double temp = sim.integrator.getTemperature();
        double Cv = energyAccumulator.getData(energyAccumulator.STANDARD_DEVIATION).getValue(0);
        Cv /= temp;
        Cv *= Cv / numAtoms;
        System.out.println("Cv " + Cv);

        // expected values based on 10^8 steps
        // for MD, avg values are very close for short and longer runs
        // stdev based on 50 x 10^6 steps with 4000 atoms (a bit larger than for 500)
        // 4 sigma should fail 1 in 16,000 runs

        double expectedP = 0.466751 + 9.74 / numAtoms;
        double stdevP = 0.005;
        if (Double.isNaN(avgP) || Math.abs(avgP - expectedP) / stdevP > 4) {
            System.exit(1);
        }

        double expectedPE = -3.81607 + 4.40 / numAtoms;
        double stdevPE = 0.0012;
        if (Double.isNaN(avgPE) || Math.abs(avgPE - expectedPE) / stdevPE > 4) {
            System.exit(2);
        }

        double expectedCv = 0.3946 + 7.724 / numAtoms;
        double stdevCv = 0.038; // stdev 500 atoms is ~2x smaller
        // at 4sigma, this isn't too useful expect that it's not super-big
        if (Double.isNaN(Cv) || Math.abs(Cv - expectedCv) / stdevCv > 4) {
            System.exit(3);
        }
    }

    public static class SimParams extends ParameterBase {
        public int numAtoms = 500;
        public int numSteps = 10000000;
        public double temperatureK = 1561.4354999999998;
        public double density = 0.08502338387498792; // V0 ==> 0 GPa
//      public double density = 0.12146197696426847; // 0.8*V0
    }
}