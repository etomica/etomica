/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.simulation.prototypes;

import etomica.action.BoxInflate;

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
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.simulation.Simulation;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;

/**
 * Simple Lennard-Jones molecular dynamics simulation in 3D
 */
public class LJMD3DNbr extends Simulation {

    public IntegratorVelocityVerlet integrator;
    public SpeciesSpheresMono species;
    public Box box;
    public P2LennardJones potential;
    public MeterPotentialEnergy energy;
    public AccumulatorAverageCollapsing avgEnergy;
    public DataPump pump;


    public LJMD3DNbr() {
        super(Space3D.getInstance());

        species = new SpeciesSpheresMono(this, space);
        species.setIsDynamic(true);
        addSpecies(species);

        PotentialMasterList potentialMaster = new PotentialMasterList(this, 4, space);
        double sigma = 1.0;
        box = this.makeBox();
        integrator = new IntegratorVelocityVerlet(this, potentialMaster, box);
        integrator.setTimeStep(0.02);

        box.setNMolecules(species, 5000);
        BoxInflate inflater = new BoxInflate(box, space);
        inflater.setTargetDensity(0.8);
        inflater.actionPerformed();

        potential = new P2LennardJones(space, sigma, 1.0);
        P2SoftSphericalTruncated p2 = new P2SoftSphericalTruncated(space, potential, 3);
        AtomType leafType = species.getLeafType();

        potentialMaster.addPotential(p2, new AtomType[]{leafType, leafType});
        integrator.getEventManager().addListener(potentialMaster.getNeighborManager(box));

        ConfigurationLattice configuration = new ConfigurationLattice(new LatticeCubicFcc(space), space);
        configuration.initializeCoordinates(box);
    }

    public static void main(String[] args) {
        final String APP_NAME = "LJMD3D";
        final LJMD3DNbr sim = new LJMD3DNbr();

sim.getController().runActivityBlocking(new etomica.action.activity.ActivityIntegrate2(sim.integrator), 300);

//        final SimulationGraphic simGraphic = new SimulationGraphic(sim, APP_NAME, 3);
//
//        ActivityIntegrate activityIntegrate = new ActivityIntegrate(sim.integrator);
//        activityIntegrate.setSleepPeriod(1);
//        sim.getController().addAction(activityIntegrate);
//
//        simGraphic.getController().getReinitButton().setPostAction(simGraphic.getPaintAction(sim.box));
//        simGraphic.getController().getDataStreamPumps().add(sim.pump);
//
//        simGraphic.makeAndDisplayFrame(APP_NAME);
//
//        DisplayTextBoxesCAE display = new DisplayTextBoxesCAE();
//        display.setAccumulator(sim.avgEnergy);
//        simGraphic.add(display);
    }

    public static class Applet extends javax.swing.JApplet {

        public void init() {
            final String APP_NAME = "LJMD3D";
            LJMD3DNbr sim = new LJMD3DNbr();
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
