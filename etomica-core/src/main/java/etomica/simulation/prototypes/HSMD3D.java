/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.simulation.prototypes;

import etomica.action.BoxImposePbc;
import etomica.action.BoxInflate;
import etomica.action.SimulationRestart;

import etomica.action.activity.ActivityIntegrate2;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DeviceNSelector;
import etomica.graphics.DisplayBox;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorHard;
import etomica.integrator.IntegratorListenerAction;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.nbr.list.NeighborListManager;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.P2HardSphere;
import etomica.potential.PotentialMaster;
import etomica.potential.PotentialMasterMonatomic;
import etomica.simulation.Simulation;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.util.ParameterBase;

/**
 * Three-dimensional hard-sphere molecular dynamics simulation, using
 * neighbor listing.
 * <p>
 * Developed as a prototype and example for the construction of a basic simulation.
 *
 * @author David Kofke and Andrew Schultz
 */
public class HSMD3D extends Simulation {

    //the following fields are made accessible for convenience to permit simple
    //mutation of the default behavior

    /**
     * The Box holding the atoms.
     */
    public final Box box;
    /**
     * The Integrator performing the dynamics.
     */
    public final IntegratorHard integrator;
    /**
     * The single hard-sphere species.
     */
    public final SpeciesSpheresMono species;
    /**
     * The hard-sphere potential governing the interactions.
     */
    public final P2HardSphere potential;

    public final PotentialMaster potentialMaster;

    /**
     * Makes a simulation using a 3D space and the default parameters.
     */
    public HSMD3D() {
        this(new HSMD3DParam());
    }

    /**
     * Makes a simulation according to the specified parameters.
     * @param params Parameters as defined by the inner class HSMD3DParam
     */
    public HSMD3D(HSMD3DParam params) {

        super(Space3D.getInstance());

        species = new SpeciesSpheresMono(this, space);
        species.setIsDynamic(true);
        addSpecies(species);

        box = this.makeBox();

        double neighborRangeFac = 1.6;
        double sigma = 1.0;
        potentialMaster = params.useNeighborLists ? new PotentialMasterList(this, sigma * neighborRangeFac, space) : new PotentialMasterMonatomic(this);

        int numAtoms = params.nAtoms;

        integrator = new IntegratorHard(this, potentialMaster, box);
        integrator.setIsothermal(false);
        integrator.setTimeStep(0.01);

        ActivityIntegrate2 ai2 = new ActivityIntegrate2(integrator);
        getController2().addActivity(ai2, Long.MAX_VALUE, 1.0);

        potential = new P2HardSphere(space, sigma, true);
        AtomType leafType = species.getLeafType();

        potentialMaster.addPotential(potential, new AtomType[]{leafType, leafType});

        box.setNMolecules(species, numAtoms);
        BoxInflate inflater = new BoxInflate(box, space);
        inflater.setTargetDensity(params.eta * 2 * space.D() / Math.PI);
        inflater.actionPerformed();
        if (space.D() == 3) {
            new ConfigurationLattice(new LatticeCubicFcc(space), space).initializeCoordinates(box);
        } else {
            new ConfigurationLattice(new LatticeOrthorhombicHexagonal(space), space).initializeCoordinates(box);
        }

        if (params.useNeighborLists) {
            NeighborListManager nbrManager = ((PotentialMasterList) potentialMaster).getNeighborManager(box);
            integrator.getEventManager().addListener(nbrManager);
        } else {
            integrator.getEventManager().addListener(new IntegratorListenerAction(new BoxImposePbc(box, space)));
        }
    }

    @Override
    public IntegratorHard getIntegrator() {
        return integrator;
    }

    /**
     * Demonstrates how this class is implemented.
     */
    public static void main(String[] args) {
        final String APP_NAME = "HSMD3D";

        HSMD3DParam params = new HSMD3DParam();
        final HSMD3D sim = new HSMD3D(params);
        final SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, APP_NAME);
        DeviceNSelector nSelector = new DeviceNSelector(sim.getController());
        nSelector.controller2 = sim.getController2();
        nSelector.setResetAction(new SimulationRestart(sim));
        nSelector.setSpecies(sim.species);
        nSelector.setBox(sim.box);

        nSelector.setPostAction(simGraphic.getPaintAction(sim.box));
        simGraphic.add(nSelector);

        simGraphic.getController().getReinitButton().setPostAction(simGraphic.getPaintAction(sim.box));

        simGraphic.makeAndDisplayFrame(APP_NAME);
        ColorSchemeByType colorScheme = ((ColorSchemeByType) ((DisplayBox) simGraphic.displayList().getFirst()).getColorScheme());
        colorScheme.setColor(sim.species.getLeafType(), java.awt.Color.red);
    }

    public static HSMD3DParam getParameters() {
        return new HSMD3DParam();
    }

    /**
     * Inner class for parameters understood by the HSMD3D constructor
     */
    public static class HSMD3DParam extends ParameterBase {
        /**
         * Number of atoms, default = 256
         */
        public int nAtoms = 256;
        /**
         * Packing fraction, default = 0.35
         */
        public double eta = 0.35;
        /**
         * Flag indicating whether neighbor list is to be used, default = true
         */
        public boolean useNeighborLists = true;
    }
}
