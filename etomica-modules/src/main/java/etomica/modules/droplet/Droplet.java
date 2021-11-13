/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.droplet;


import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.potential.BondingInfo;
import etomica.potential.PotentialMaster;
import etomica.potential.compute.PotentialComputeAggregate;
import etomica.potential.compute.PotentialComputeField;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularNonperiodic;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesGeneral;

/**
 * Mesoscale simulation for Droplet module.
 *
 * @author Andrew Schultz
 */
public class Droplet extends Simulation {

    public final SpeciesGeneral species;
    public final Box box;
    public final IntegratorDroplet integrator;
    public final P2Cohesion p2;
    public final P1Smash p1Smash;
    public final ConfigurationDroplet config;
    public final AtomTestLiquid liquidFilter;
    public final MeterDeformation meterDeformation;

    public Droplet() {
        super(Space3D.getInstance());

        //species
        species = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this), true);
        addSpecies(species);

        box = this.makeBox(new BoundaryRectangularNonperiodic(space));
        int numAtoms = 2000;
        PotentialMaster potentialMaster = new PotentialMaster(getSpeciesManager(), box, BondingInfo.noBonding());
        PotentialComputeField pcField = new PotentialComputeField(getSpeciesManager(), box);
        PotentialComputeAggregate pcAgg = new PotentialComputeAggregate(potentialMaster, pcField);

        //controller and integrator
        integrator = new IntegratorDroplet(this, pcAgg, box);
        getController().addActivity(new ActivityIntegrate(integrator), Long.MAX_VALUE, 0.0);
        integrator.setTimeStep(0.2);
        integrator.setTemperature(0);

        //potentials
        AtomType leafType = species.getLeafType();

        p2 = new P2Cohesion(space);
        p2.setEpsilon(1.0);
        double vol = 4.0 / 3.0 * Math.PI;
        p2.setDv(vol / numAtoms);
        potentialMaster.setPairPotential(leafType, leafType, p2);

        p1Smash = new P1Smash(space);
        p1Smash.setG(1.0);
        pcField.setFieldPotential(leafType, p1Smash);

        //construct box
        Vector dim = space.makeVector();
        dim.E(new double[]{4, 4, 4});
        box.getBoundary().setBoxSize(dim);
        box.setNMolecules(species, numAtoms);

        config = new ConfigurationDroplet(random, space);
        config.initializeCoordinates(box);

        meterDeformation = new MeterDeformation(space);
        meterDeformation.setBox(box);

        liquidFilter = new AtomTestLiquid(space, meterDeformation);
        liquidFilter.setCutoff(0.9);

        p2.setLiquidFilter(liquidFilter);
    }

    @Override
    public IntegratorDroplet getIntegrator() {
        return integrator;
    }

    public static void main(String[] args) {
        Droplet sim = new Droplet();
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, Long.MAX_VALUE));
    }//end of main
}
