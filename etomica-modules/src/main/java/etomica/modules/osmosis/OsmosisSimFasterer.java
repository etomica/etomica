/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.osmosis;

import etomica.action.IAction;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorHardFasterer;
import etomica.integrator.IntegratorMDFasterer;
import etomica.lattice.LatticeCubicSimple;
import etomica.math.geometry.Plane;
import etomica.potential.P1HardBoundary;
import etomica.potential.P2HardGeneric;
import etomica.potential.P2HardSphere;
import etomica.potential.compute.NeighborManagerSimpleHard;
import etomica.potential.compute.PotentialComputeField;
import etomica.potential.compute.PotentialComputePair;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularNonperiodic;
import etomica.space.Space;
import etomica.space2d.Space2D;
import etomica.space2d.Vector2D;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.species.SpeciesGeneral;

import java.awt.*;

/**
 * Osmosis simulation.
 *
 * @author Jhumpa Adhikari
 * @author Andrew Schultz
 */

public class OsmosisSimFasterer extends Simulation {

    protected final static int initialSolute = 10;
    protected final static int initialSolvent = 50;
    public IntegratorHardFasterer integrator;
    public SpeciesGeneral speciesSolvent, speciesSolute;
    public Box box;
    public P1HardBoundary boundaryHardA;
    public P1HardBoundary boundaryHardB;


    public OsmosisSimFasterer(Space _space) {
        super(_space);

        final double sigma = 1.0;

        speciesSolvent = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this), true);
        addSpecies(speciesSolvent);
        speciesSolute = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this), true);
        addSpecies(speciesSolute);

        //construct box
        box = this.makeBox(new BoundaryRectangularNonperiodic(space));

        NeighborManagerSimpleHard neighborManager = new NeighborManagerSimpleHard(box);
        PotentialComputePair potentialMaster = new PotentialComputePair(getSpeciesManager(), box, neighborManager);
        P2HardGeneric potential = P2HardSphere.makePotential(sigma);
        potentialMaster.setPairPotential(speciesSolvent.getLeafType(), speciesSolvent.getLeafType(), potential);
        potentialMaster.setPairPotential(speciesSolute.getLeafType(), speciesSolvent.getLeafType(), potential);
        potentialMaster.setPairPotential(speciesSolute.getLeafType(), speciesSolute.getLeafType(), potential);

        PotentialComputeField pcField = new PotentialComputeField(getSpeciesManager(), box);
        //Boundary potential for the solvent
        boundaryHardA = new P1HardBoundary(space, false, box);
        pcField.setFieldPotential(speciesSolvent.getLeafType(), boundaryHardA);
        boundaryHardA.setCollisionRadius(0.5 * sigma);

        //Boundary potential for the solute
        boundaryHardB = new P1HardBoundaryOsmosis(box);
        pcField.setFieldPotential(speciesSolute.getLeafType(), boundaryHardB);
        boundaryHardB.setCollisionRadius(0.5 * sigma);

        if (space instanceof Space2D) { // 2D
            box.getBoundary().setBoxSize(new Vector2D(10.0, 10.0));
        } else if (space instanceof Space3D) { // 3D
            box.getBoundary().setBoxSize(new Vector3D(10.0, 10.0, 10.0));
        }
        box.setNMolecules(speciesSolvent, initialSolvent);
        box.setNMolecules(speciesSolute, initialSolute);

        integrator = new IntegratorHardFasterer(IntegratorHardFasterer.extractHardPotentials(potentialMaster), IntegratorHardFasterer.extractFieldPotentials(pcField),
                neighborManager, random, 0.05, 1.0, box, getSpeciesManager(), null);
        integrator.setThermostat(IntegratorMDFasterer.ThermostatType.ANDERSEN_SINGLE);

        ConfigurationLattice config = new ConfigurationLattice(new LatticeCubicSimple(space, 1.0), space);
        config.initializeCoordinates(box);
    }


    public static void main(String[] args) {
        Space sp = Space3D.getInstance();
        final OsmosisSimFasterer sim = new OsmosisSimFasterer(sp);
        final ConfigurationLattice config = new ConfigurationLattice(new LatticeCubicSimple(sp, 1.0), sp);
        config.initializeCoordinates(sim.box);
        Plane plane = new Plane(sim.getSpace());

        final SimulationGraphic simGraphic = new SimulationGraphic(sim, "Osmosis Sim");
        ((etomica.graphics.DisplayBoxCanvasG3DSys) simGraphic.getDisplayBox(sim.box).canvas).addPlane(plane);
        simGraphic.getController().getReinitButton().setPostAction(new IAction() {
            public void actionPerformed() {
                config.initializeCoordinates(sim.box);
                simGraphic.getDisplayBox(sim.box).repaint();
            }
        });
        simGraphic.makeAndDisplayFrame("Osmosis Sim");
        ColorSchemeByType colorScheme = new ColorSchemeByType();
        colorScheme.setColor(sim.speciesSolvent.getLeafType(), Color.blue);
        colorScheme.setColor(sim.speciesSolute.getLeafType(), Color.red);
        simGraphic.getDisplayBox(sim.box).setColorScheme(colorScheme);
        config.initializeCoordinates(sim.box);
        simGraphic.getDisplayBox(sim.box).repaint();
        sim.integrator.setTimeStep(0.05);
        sim.getController().setSleepPeriod(0);
        sim.getController().addActivity(new ActivityIntegrate(sim.integrator));
    }

} 
