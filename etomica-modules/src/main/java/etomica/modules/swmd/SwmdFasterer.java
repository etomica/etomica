/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.swmd;

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
import etomica.potential.compute.PotentialComputePairGeneral;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesGeneral;
import etomica.units.*;

public class SwmdFasterer extends Simulation {

    public SpeciesGeneral species;
    public Box box;
    public IntegratorHardFasterer integrator;
    public P2HardGeneric p2sqw;

    public SwmdFasterer(Space _space) {
        super(_space);

        //species
        species = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this), true);//index 1
        ((ElementSimple) species.getLeafType().getElement()).setMass(Dalton.UNIT.toSim(space.D() == 3 ? 131 : 40));
        addSpecies(species);
        box = this.makeBox();
        NeighborManagerSimpleHard neighborManager = new NeighborManagerSimpleHard(box);

        PotentialComputePairGeneral potentialMaster = new PotentialComputePairGeneral(getSpeciesManager(), box, neighborManager);

        int N = space.D() == 3 ? 256 : 100;  //number of atoms

        double sigma = 4.0;
        double lambda = 2.0;
        double epsilon = new UnitRatio(Joule.UNIT, Mole.UNIT).toSim(space.D() == 3 ? 1000 : 1500);
        p2sqw = new P2HardGeneric(new double[]{sigma, sigma * lambda}, new double[]{Double.POSITIVE_INFINITY, -epsilon}, true);
        potentialMaster.setPairPotential(species.getLeafType(), species.getLeafType(), p2sqw);

        //controller and integrator
        integrator = new IntegratorHardFasterer(IntegratorHardFasterer.extractHardPotentials(potentialMaster), neighborManager,
                random, 1.0, Kelvin.UNIT.toSim(300), box);
        integrator.setIsothermal(false);
        integrator.setThermostat(IntegratorMDFasterer.ThermostatType.ANDERSEN_SINGLE);
        integrator.setThermostatInterval(1);
        getController().addActivity(new ActivityIntegrate(integrator));

        integrator.setMaxCollisionDiameter(species.getLeafType(), sigma * lambda);

        //construct box
        Vector dim = space.makeVector();
        dim.E(space.D() == 3 ? 30 : 50);
        box.getBoundary().setBoxSize(dim);
        box.setNMolecules(species, N);
        new ConfigurationLattice(space.D() == 3 ? new LatticeCubicFcc(space) : new LatticeOrthorhombicHexagonal(space), space).initializeCoordinates(box);

        integrator.getEventManager().addListener(new IntegratorListenerAction(new BoxImposePbc(box, space)));
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

        SwmdFasterer sim = new SwmdFasterer(space);
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, Long.MAX_VALUE));
    }
}
