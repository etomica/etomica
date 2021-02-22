/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.insertion;

import etomica.action.BoxImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.chem.elements.ElementSimple;
import etomica.config.ConfigurationLattice;
import etomica.integrator.IntegratorHardFasterer;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorMDFasterer;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.potential.P2HardGeneric;
import etomica.potential.compute.NeighborManagerSimpleHard;
import etomica.potential.compute.PotentialComputePair;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesGeneral;
import etomica.util.random.RandomMersenneTwister;

public class InsertionFasterer extends Simulation {

    public final SpeciesGeneral species, speciesGhost;
    public final Box box;
    public final IntegratorHardFasterer integrator;
    public final P2HardGeneric p2sqw;
    public final P2HardGeneric potentialGhost;


    public InsertionFasterer(Space _space) {
        super(_space);
        setRandom(new RandomMersenneTwister(2));

        // species
        species = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this), true);//index 1
        addSpecies(species);
        speciesGhost = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this), true);
        ((ElementSimple) speciesGhost.getLeafType().getElement()).setMass(Double.POSITIVE_INFINITY);
        addSpecies(speciesGhost);

        int N = space.D() == 3 ? 256 : 100;  //number of atoms

        double sigma = 1.0;
        double lambda = 1.5;

        //controller and integrator
        box = this.makeBox();
        NeighborManagerSimpleHard neighborManager = new NeighborManagerSimpleHard(box);
        PotentialComputePair potentialCompute = new PotentialComputePair(this, box, neighborManager);

        //instantiate several potentials for selection in combo-box
        p2sqw = new P2HardGeneric(new double[]{sigma, sigma * lambda}, new double[]{Double.POSITIVE_INFINITY, -1.0}, true);
        potentialCompute.setPairPotential(species.getLeafType(), species.getLeafType(), p2sqw);

        potentialGhost = new P2HardGeneric(new double[]{sigma, sigma * lambda}, new double[]{0, 0});
        potentialCompute.setPairPotential(species.getLeafType(), speciesGhost.getLeafType(), potentialGhost);

        integrator = new IntegratorHardFasterer(IntegratorHardFasterer.extractHardPotentials(potentialCompute), neighborManager, random, 0.2, 1.0, box);
        integrator.setIsothermal(false);
        integrator.setThermostat(IntegratorMDFasterer.ThermostatType.ANDERSEN_SCALING);
        integrator.setThermostatNoDrift(true);
        integrator.setThermostatInterval(1);

        //construct box
        Vector dim = space.makeVector();
        dim.E(space.D() == 3 ? 8 : 13.5);
        box.getBoundary().setBoxSize(dim);
        box.setNMolecules(species, N);
        new ConfigurationLattice(space.D() == 3 ? new LatticeCubicFcc(space) : new LatticeOrthorhombicHexagonal(space), space).initializeCoordinates(box);
        box.setNMolecules(speciesGhost, 1);

        integrator.getEventManager().addListener(new IntegratorListenerAction(new BoxImposePbc(box, space)));
        potentialCompute.init();
    }

    public IntegratorHardFasterer getIntegrator() {
        return integrator;
    }

    public static void main(String[] args) {
        Space space = Space3D.getInstance();
        if (args.length != 0) {
            try {
                int D = Integer.parseInt(args[0]);
                if (D == 3) {
                    space = Space3D.getInstance();
                }
            } catch (NumberFormatException e) {
            }
        }

        InsertionFasterer sim = new InsertionFasterer(space);
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, Long.MAX_VALUE));
    }
}
