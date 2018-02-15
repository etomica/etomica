/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.simulation.prototypes;

import etomica.action.BoxInflate;
import etomica.action.activity.ActivityIntegrate;
import etomica.action.activity.Controller;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorAverageCollapsing;
import etomica.data.DataPump;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.graphics.DisplayTextBoxesCAE;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.lattice.LatticeCubicFcc;
import etomica.nbr.list.PotentialMasterListFast;
import etomica.nbr.list.PotentialMasterListFaster;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.potential.PotentialCalculationForceSum;
import etomica.simulation.Simulation;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;

/**
 * Simple Lennard-Jones molecular dynamics simulation in 3D
 */
public class LJMD3DFaster extends Simulation {

    public IntegratorVelocityVerlet integrator;
    public SpeciesSpheresMono species;
    public Box box;
    public P2LennardJones potential;
    public Controller controller;
    public MeterPotentialEnergy energy;
    public AccumulatorAverageCollapsing avgEnergy;
    public DataPump pump;

    public LJMD3DFaster() {
        super(Space3D.getInstance());
        PotentialMasterListFaster potentialMaster = new PotentialMasterListFaster(this, 4, space);
        potentialMaster.lrcMaster().setEnabled(false);
        double sigma = 1.0;
        integrator = new IntegratorVelocityVerlet(this, potentialMaster, space);
        integrator.setTimeStep(0.02);
        integrator.setForceSum(new PotentialCalculationForceSum());
        ActivityIntegrate activityIntegrate = new ActivityIntegrate(integrator);
        activityIntegrate.setSleepPeriod(1);
        getController().addAction(activityIntegrate);
        species = new SpeciesSpheresMono(this, space);
        species.setIsDynamic(true);
        addSpecies(species);
        box = new Box(space);
        addBox(box);
        box.setNMolecules(species, 5000);
        BoxInflate inflater = new BoxInflate(box, space);
        inflater.setTargetDensity(0.8);
        inflater.actionPerformed();
        potential = new P2LennardJones(space, sigma, 1.0);
        P2SoftSphericalTruncated p2 = new P2SoftSphericalTruncated(space, potential, 3);
        AtomType leafType = species.getLeafType();

        potentialMaster.addPotential(p2, new AtomType[]{leafType, leafType});
        integrator.setBox(box);
        integrator.getEventManager().addListener(potentialMaster.getNeighborManager(box));


        ConfigurationLattice configuration = new ConfigurationLattice(new LatticeCubicFcc(space), space);
        configuration.initializeCoordinates(box);
    }

    public static void main(String[] args) {
        final String APP_NAME = "LJMD3D";
        final LJMD3DFaster sim = new LJMD3DFaster();
        final SimulationGraphic simGraphic = new SimulationGraphic(sim, APP_NAME, 3);

        simGraphic.getController().getReinitButton().setPostAction(simGraphic.getPaintAction(sim.box));
        simGraphic.getController().getDataStreamPumps().add(sim.pump);

        simGraphic.makeAndDisplayFrame(APP_NAME);

        DisplayTextBoxesCAE display = new DisplayTextBoxesCAE();
        display.setAccumulator(sim.avgEnergy);
        simGraphic.add(display);
    }

    public static class Applet extends javax.swing.JApplet {

        public void init() {
            final String APP_NAME = "LJMD3D";
            LJMD3DFaster sim = new LJMD3DFaster();
            final SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.GRAPHIC_ONLY, APP_NAME, 3);

            simGraphic.getController().getReinitButton().setPostAction(simGraphic.getPaintAction(sim.box));
            simGraphic.getController().getDataStreamPumps().add(sim.pump);

            DisplayTextBoxesCAE display = new DisplayTextBoxesCAE();
            display.setAccumulator(sim.avgEnergy);
            simGraphic.add(display);
            getContentPane().add(simGraphic.getPanel());
        }
    }

}
