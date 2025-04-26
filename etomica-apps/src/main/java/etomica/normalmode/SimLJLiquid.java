/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.action.BoxInflate;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomSourceFixed0;
import etomica.atom.AtomType;
import etomica.atom.DiameterHashByType;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.DataPumpListener;
import etomica.data.DataSourceCountSteps;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.meter.MeterPressure;
import etomica.data.meter.MeterPressureFromIntegrator;
import etomica.graphics.*;
import etomica.integrator.*;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.lattice.LatticeCubicBcc;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.LatticeCubicSimple;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.nbr.list.NeighborListManager;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.*;
import etomica.potential.compute.PotentialCompute;
import etomica.potential.compute.PotentialComputeAggregate;
import etomica.potential.compute.PotentialComputeField;
import etomica.potential.compute.PotentialComputePair;
import etomica.simulation.Simulation;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesGeneral;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

import java.awt.*;

/**
 * Simple Lennard-Jones Monte Carlo simulation in 3D.
 */
public class SimLJLiquid extends Simulation {
    public IntegratorMD integrator;
    public PotentialComputePair pm;
    public SpeciesGeneral species;
    public P2LennardJones potential;
    public double temperature;
    public Box box;

    public SimLJLiquid(int nAtoms, double rho, double temperature, double rc, boolean isLRC, double dt, double gamma, double nu, Thermostat thermostat) {
        super(Space3D.getInstance());
        species = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this), true);
        AtomType leafType = species.getLeafType();
        addSpecies(species);
        box = this.makeBox();
        box.setNMolecules(species, nAtoms);
        BoxInflate inflater = new BoxInflate(box, space);
        inflater.setTargetDensity(rho);
        inflater.actionPerformed();

        ConfigurationLattice config = new ConfigurationLattice(new LatticeCubicSimple(space), space);//BCC
        config.initializeCoordinates(box);

        double sigma = 1, eps = 1;
        potential = new P2LennardJones(sigma, eps);
        P2SoftSphericalTruncatedForceShifted potentialTruncated = new P2SoftSphericalTruncatedForceShifted(potential, rc);
        PotentialMasterList pm = new PotentialMasterList(this.getSpeciesManager(), box, 2, 1.5*rc, BondingInfo.noBonding());
        pm.setPairPotential(leafType, leafType, potentialTruncated);

        if (!isLRC) {
            pm.doAllTruncationCorrection = false;
        }

        if (thermostat == Thermostat.AT) { //Andersen Thermostat
            double cdf = 1.0 - Math.exp(-nu*dt);
            integrator = new IntegratorVelocityVerlet(pm, random, dt, temperature, box);
            integrator.setThermostatInterval(1);
            integrator.setCdf_andersen(cdf);
            System.out.println(" AT thermostat (nu: " + nu + " , cdf: " + cdf + ")");
        } else if (thermostat == Thermostat.LT) { //Langevin Thermostat
            integrator = new IntegratorLangevin(pm, random, dt, temperature, box, gamma);
            integrator.setThermostatNoDrift(false);
            System.out.println(" LT thermostat (gamma: " + gamma+")");
        } else if (thermostat == Thermostat.BD) { //Brownian Sampling
            integrator = new IntegratorBrownian(pm, random, dt, temperature, box);
            integrator.setThermostatNoDrift(false);
            System.out.println(" BD thermostat (h: " + dt+")");
        }
        integrator.setIsothermal(true);
    }

    public static void main(String[] args) {
        long t1 = System.nanoTime();
        SimParams params = new SimParams();
        ParseArgs.doParseArgs(params, args);
        boolean isGraphic = params.isGraphic;
        int nAtoms = params.nAtoms;
        int nSteps = params.nSteps;
        int nStepsEq = nSteps/10;
        double temperature = params.temperature;
        double rho = params.rho;
        double rc = params.rc;
        double dt = params.dt;
        Thermostat thermostat = params.thermostat;
        boolean isLRC = params.isLRC;
        double gamma = params.gamma;
        double nu = params.nu;

        SimLJLiquid sim = new SimLJLiquid(nAtoms, rho, temperature, rc, isLRC, dt, gamma, nu, thermostat);
        double L = sim.box.getBoundary().getBoxSize().getX(0);
        System.out.println(" isLRC: " + isLRC);
        System.out.println(" N: " + nAtoms + " , density: " + rho + " , L/2: " + L/2 + " , T: " + temperature);
        System.out.println(" rc: " + rc);
        System.out.println(" dt: " + dt);

        if (isGraphic) {
            sim.getController().addActivity(new ActivityIntegrate(sim.integrator));
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
            simGraphic.setPaintInterval(sim.box, 1);
            ColorScheme colorScheme = new ColorScheme() {
                protected Color[] allColors;

                public Color getAtomColor(IAtom a) {
                    if (allColors==null) {
                        allColors = new Color[768];
                        for (int i=0; i<256; i++) {
                            allColors[i] = new Color(255-i,i,0);
                        }
                        for (int i=0; i<256; i++) {
                            allColors[i+256] = new Color(0,255-i,i);
                        }
                        for (int i=0; i<256; i++) {
                            allColors[i+512] = new Color(i,0,255-i);
                        }
                    }
                    return allColors[(2*a.getLeafIndex()) % 768];
                }
            };
            simGraphic.getDisplayBox(sim.box).setColorScheme(colorScheme);

            ((DiameterHashByType)((DisplayBox)simGraphic.displayList().getFirst()).getDiameterHash()).setDiameter(sim.species().getAtomType(0),0.9);

            DisplayTextBox timer = new DisplayTextBox();
            DataSourceCountSteps counter = new DataSourceCountSteps(sim.integrator);
            DataPumpListener counterPump = new DataPumpListener(counter, timer, 100);
            sim.integrator.getEventManager().addListener(counterPump);
            simGraphic.getPanel().controlPanel.add(timer.graphic());

            DisplayPlot pePlot = new DisplayPlot();
            simGraphic.add(pePlot);

            simGraphic.makeAndDisplayFrame();
            return;
        }
        System.out.flush();

        int intervalSample = 10;
        System.out.println(" intervalSample: " + intervalSample);
        int blocks = 100;
        long blockSize = params.nSteps / (intervalSample * blocks);
        System.out.println(" steps: " +  nSteps);
        System.out.println(" nBlocks: " + blocks + " blocksize: " + blockSize + " interval: " + intervalSample);
        System.out.println();

        // Equilibaration
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, nStepsEq));
        System.out.println(" Done with " + nStepsEq + " steps equilibration.");

        MeterPotentialEnergyFromIntegrator meterU = new MeterPotentialEnergyFromIntegrator(sim.integrator);
        AccumulatorAverageFixed accumulatorU = new AccumulatorAverageFixed(blockSize);
        DataPumpListener accumulatorPumpU = new DataPumpListener(meterU, accumulatorU, intervalSample);
        sim.integrator.getEventManager().addListener(accumulatorPumpU);

        MeterPressureFromIntegrator meterP = new MeterPressureFromIntegrator(sim.integrator);
        AccumulatorAverage accumulatorP = new AccumulatorAverageFixed(blockSize);
        DataPumpListener pPump = new DataPumpListener(meterP, accumulatorP, intervalSample);
        sim.integrator.getEventManager().addListener(pPump);

        // Production
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, nSteps));

        double avgU = accumulatorU.getData(accumulatorU.AVERAGE).getValue(0);
        double errU = accumulatorU.getData(accumulatorU.ERROR).getValue(0);
        double corU = accumulatorU.getData(accumulatorU.BLOCK_CORRELATION).getValue(0);
        System.out.println();
        System.out.println(" U: " + avgU/nAtoms + "  " + errU/nAtoms + "    " + corU);

        double avgP = accumulatorP.getData(accumulatorP.AVERAGE).getValue(0);
        double errP = accumulatorP.getData(accumulatorP.ERROR).getValue(0);
        double corP = accumulatorP.getData(accumulatorP.BLOCK_CORRELATION).getValue(0);
        System.out.println(" P: " + avgP + "  " + errP + "    " + corP);

        long t2 = System.nanoTime();
        System.out.println("\n time: " + (t2 - t1)/1.0e9/60.0 + " min");
    }

    public enum Thermostat {LT, AT, BD};

    public static class SimParams extends ParameterBase {
        public boolean isGraphic = !false;
        public boolean isLRC = true;
        public int nAtoms = 500;
        public int nSteps = 100000;
        public double temperature = 1.0;
        public double rho = 0.8;
        public double rc = 2.5;
        public double dt = 0.0005;
        public Thermostat thermostat = Thermostat.BD;
        public double gamma = 200; //LT
        public double nu = 50; // AT
    }
}