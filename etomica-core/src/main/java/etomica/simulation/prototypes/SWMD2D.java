/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.simulation.prototypes;

import etomica.action.BoxImposePbc;

import etomica.action.activity.ActivityIntegrate2;
import etomica.action.activity.Controller;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.graphics.DisplayBox;
import etomica.integrator.IntegratorHard;
import etomica.integrator.IntegratorListenerAction;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.potential.P2SquareWell;
import etomica.potential.PotentialMaster;
import etomica.potential.PotentialMasterMonatomic;
import etomica.simulation.Simulation;
import etomica.space2d.Space2D;
import etomica.species.SpeciesSpheresMono;

/**
 * Simple square-well molecular dynamics simulation in 2D
 */
public class SWMD2D extends Simulation {

    private static final long serialVersionUID = 1L;
    public IntegratorHard integrator;
    public SpeciesSpheresMono species;
    public Box box;
    public P2SquareWell potential;
    public DisplayBox display;

    public SWMD2D() {
        super(Space2D.getInstance());

        species = new SpeciesSpheresMono(this, space);
        species.setIsDynamic(true);
        addSpecies(species);

        PotentialMaster potentialMaster = new PotentialMasterMonatomic(this);
        double sigma = 0.8;
        box = this.makeBox();
        integrator = new IntegratorHard(this, potentialMaster, box);
        integrator.setTimeStep(0.02);
        integrator.setTemperature(1.);
        getController2().addActivity(new ActivityIntegrate2(integrator)).setSleepPeriod(1);
        AtomType leafType = species.getLeafType();
        box.setNMolecules(species, 50);
        new ConfigurationLattice(new LatticeOrthorhombicHexagonal(space), space).initializeCoordinates(box);
        potential = new P2SquareWell(space);
        potential.setCoreDiameter(sigma);
        potentialMaster.addPotential(potential, new AtomType[]{leafType, leafType});
        integrator.getEventManager().addListener(new IntegratorListenerAction(new BoxImposePbc(box, space)));
    }
}
