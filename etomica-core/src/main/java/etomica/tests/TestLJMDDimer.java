/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.tests;

import etomica.action.ActionIntegrate;
import etomica.action.BoxImposePbc;
import etomica.action.BoxInflate;
import etomica.action.activity.ActivityIntegrate;
import etomica.action.activity.Controller;
import etomica.atom.AtomType;
import etomica.atom.iterator.ApiIndexList;
import etomica.atom.iterator.AtomIteratorBasisDependent;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.config.ConformationLinear;
import etomica.data.AccumulatorAverageCollapsing;
import etomica.data.DataPump;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.graphics.ColorSchemeRandomByMolecule;
import etomica.graphics.DisplayTextBoxesCAE;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.lattice.LatticeCubicFcc;
import etomica.potential.*;
import etomica.simulation.Simulation;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesSpheres;
import etomica.species.SpeciesSpheresMono;

import java.util.zip.Inflater;

/**
 * Simple Lennard-Jones molecular dynamics simulation in 3D
 */
public class TestLJMDDimer extends Simulation {

    public IntegratorVelocityVerlet integrator;
    public SpeciesSpheres species;
    public Box box;
    public P2LennardJones potential;
    public Controller controller;
    public MeterPotentialEnergy energy;
    public AccumulatorAverageCollapsing avgEnergy;
    public DataPump pump;


    public TestLJMDDimer() {
        super(Space3D.getInstance());

        species = new SpeciesSpheres(this, space, 2);
        species.setConformation(new ConformationLinear(space, 0.5));
        species.setIsDynamic(true);
        addSpecies(species);

        PotentialMaster potentialMaster = new PotentialMaster();
        double sigma = 1.0;
        box = this.makeBox();
        integrator = new IntegratorVelocityVerlet(this, potentialMaster, box);
        integrator.setTimeStep(0.005);
        integrator.setTemperature(2);
        integrator.setIsothermal(true);
        ActivityIntegrate activityIntegrate = new ActivityIntegrate(integrator);
        activityIntegrate.setSleepPeriod(0);
        activityIntegrate.setMaxSteps(1000);
        getController().addAction(activityIntegrate);
        box.setNMolecules(species, 512);
        new BoxInflate(box, space, 0.5).actionPerformed();
        System.out.println("box size: "+box.getBoundary().getBoxSize());

        potential = new P2LennardJones(space, sigma, 1.0);
        AtomType leafType = species.getLeafType();
        P2SoftSphericalTruncatedForceShifted p2 = new P2SoftSphericalTruncatedForceShifted(space, potential, 3.0);
        potentialMaster.addPotential(p2, new AtomType[]{leafType, leafType});

        P2Harmonic pBond = new P2Harmonic(space, 100, 0.51);
        PotentialGroup p1 = new PotentialGroup(1, space);
        ApiIndexList bondIterator = new ApiIndexList(new int[][]{{0,1}});
        p1.addPotential(pBond, bondIterator);
        potentialMaster.addPotential(p1, new ISpecies[]{species});

        BoxImposePbc imposepbc = new BoxImposePbc(space);
        imposepbc.setBox(box);
        integrator.getEventManager().addListener(new IntegratorListenerAction(imposepbc));

        ConfigurationLattice configuration = new ConfigurationLattice(new LatticeCubicFcc(space), space);
        configuration.initializeCoordinates(box);
        energy = new MeterPotentialEnergy(potentialMaster, box);
        System.out.println("u0: "+energy.getDataAsScalar());
        avgEnergy = new AccumulatorAverageCollapsing();
        avgEnergy.setPushInterval(10);
        pump = new DataPump(energy, avgEnergy);
        IntegratorListenerAction pumpListener = new IntegratorListenerAction(pump);
        pumpListener.setInterval(10);
        integrator.getEventManager().addListener(pumpListener);
    }

    public static void main(String[] args) {
        final String APP_NAME = "LJMDDimer";
        final TestLJMDDimer sim = new TestLJMDDimer();
        long t0 = System.nanoTime();
        sim.getController().actionPerformed();
        long t1 = System.nanoTime();
        System.out.println((t1 - t0) / 1e6);
//        final SimulationGraphic simGraphic = new SimulationGraphic(sim, APP_NAME, 3);
//
//        simGraphic.getController().getReinitButton().setPostAction(simGraphic.getPaintAction(sim.box));
//        simGraphic.getController().getDataStreamPumps().add(sim.pump);
//        simGraphic.getDisplayBox(sim.box).setColorScheme(new ColorSchemeRandomByMolecule(sim, sim.box, sim.getRandom()));
//
//        simGraphic.makeAndDisplayFrame(APP_NAME);
//
//        DisplayTextBoxesCAE display = new DisplayTextBoxesCAE();
//        display.setAccumulator(sim.avgEnergy);
//        simGraphic.add(display);
    }
}
