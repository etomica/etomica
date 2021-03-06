/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.simulation.prototypes;


import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.data.meter.MeterEnergy;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayPlot;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.potential.P2LennardJones;
import etomica.potential.PotentialMaster;
import etomica.potential.PotentialMasterMonatomic;
import etomica.simulation.Simulation;
import etomica.space2d.Space2D;
import etomica.species.SpeciesGeneral;

/**
 * Simple Lennard-Jones molecular dynamics simulation in 2D
 */
public class LJMD2D extends Simulation {

    private static final long serialVersionUID = 1L;
    public IntegratorVelocityVerlet integrator;
    public SpeciesGeneral species;
    public Box box;
    public P2LennardJones potential;
    public DisplayBox display;
    public DisplayPlot plot;
    public MeterEnergy energy;

    public LJMD2D() {
        super(Space2D.getInstance());

        species = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this), true);
        addSpecies(species);

        PotentialMaster potentialMaster = new PotentialMasterMonatomic(this);
        box = this.makeBox();
        integrator = new IntegratorVelocityVerlet(this, potentialMaster, box);
        integrator.setTimeStep(0.01);
        getController().setSleepPeriod(2);
        getController().addActivity(new ActivityIntegrate(integrator));
        box.setNMolecules(species, 50);
        new ConfigurationLattice(new LatticeOrthorhombicHexagonal(space), space).initializeCoordinates(box);
        potential = new P2LennardJones(space);
        potentialMaster.addPotential(potential, new AtomType[]{species.getLeafType(), species.getLeafType()});

        energy = new MeterEnergy(potentialMaster, box);

    }

}
