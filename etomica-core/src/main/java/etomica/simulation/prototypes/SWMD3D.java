/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

//Source file generated by Etomica

package etomica.simulation.prototypes;

import etomica.action.BoxInflate;
import etomica.action.activity.ActivityIntegrate;
import etomica.action.activity.Controller;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.data.DataPumpListener;
import etomica.data.meter.MeterPressureHard;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayTextBox;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorHard;
import etomica.lattice.LatticeCubicFcc;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.P2SquareWell;
import etomica.simulation.Simulation;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;

/**
 * Simple square-well molecular dynamics simulation in 3D
 */
public class SWMD3D extends Simulation {

    public IntegratorHard integrator;
    public SpeciesSpheresMono species;
    public Box box;
    public P2SquareWell potential;
    public Controller controller;
    public DisplayBox display;

    public SWMD3D() {
        super(Space3D.getInstance());
        PotentialMasterList potentialMaster = new PotentialMasterList(this, 2.5, space);

        box = new Box(space);
        integrator = new IntegratorHard(this, potentialMaster, space, box);
        integrator.setTimeStep(0.01);
        integrator.setIsothermal(true);
        integrator.setTemperature(1);
        double lambda = 2;
        ActivityIntegrate activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);
        addBox(box);
        potential = new P2SquareWell(space);
        potential.setLambda(lambda);

        species = new SpeciesSpheresMono(this, space);
        species.setIsDynamic(true);
        addSpecies(species);
        box.setNMolecules(species, 108);

        potentialMaster.addPotential(potential, new AtomType[]{species.getLeafType(), species.getLeafType()});

        integrator.setBox(box);
        integrator.getEventManager().addListener(potentialMaster.getNeighborManager(box));

        BoxInflate inflater = new BoxInflate(box, space);
        inflater.setTargetDensity(0.0405);
        inflater.actionPerformed();
        ConfigurationLattice configuration = new ConfigurationLattice(new LatticeCubicFcc(space), space);
        configuration.initializeCoordinates(box);
    }

    /**
     * Demonstrates how this class is implemented.
     */
    public static void main(String[] args) {
        final String APP_NAME = "SWMD3D";

        final SWMD3D sim = new SWMD3D();
        final SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, APP_NAME);

        simGraphic.getController().getReinitButton().setPostAction(simGraphic.getPaintAction(sim.box));

        simGraphic.makeAndDisplayFrame(APP_NAME);
        ColorSchemeByType colorScheme = ((ColorSchemeByType) ((DisplayBox) simGraphic.displayList().getFirst()).getColorScheme());
        colorScheme.setColor(sim.species.getLeafType(), java.awt.Color.red);

        MeterPressureHard pMeter = new MeterPressureHard(sim.space);
        pMeter.setIntegrator(sim.integrator);

        DisplayTextBox pdisplay = new DisplayTextBox();
        DataPumpListener pPump = new DataPumpListener(pMeter, pdisplay, 100);
        sim.integrator.getEventManager().addListener(pPump);
        simGraphic.add(pdisplay);
    }
}
