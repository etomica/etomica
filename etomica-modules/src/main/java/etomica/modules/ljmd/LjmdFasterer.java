/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.ljmd;


import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.integrator.IntegratorMDFasterer;
import etomica.integrator.IntegratorVelocityVerletFasterer;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.nbr.list.PotentialMasterListFasterer;
import etomica.potential.BondingInfo;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space2d.Space2D;
import etomica.space3d.Space3D;
import etomica.species.SpeciesGeneral;

public class LjmdFasterer extends Simulation {

    public SpeciesGeneral species;
    public Box box;
    public IntegratorVelocityVerletFasterer integrator;


    public LjmdFasterer(Space _space) {
        super(_space);

        //species
        species = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this), true);//index 1
        addSpecies(species);

        box = this.makeBox();
        PotentialMasterListFasterer potentialMaster = new PotentialMasterListFasterer(getSpeciesManager(), box, 2, 2.99, BondingInfo.noBonding());

        int N = 182;  //number of atoms

        //controller and integrator
        integrator = new IntegratorVelocityVerletFasterer(potentialMaster, random, 0.01, 1.0, box);
        integrator.setIsothermal(false);
        integrator.setThermostat(IntegratorMDFasterer.ThermostatType.ANDERSEN_SINGLE);
        integrator.setThermostatNoDrift(true);
        integrator.setThermostatInterval(1);
        //   integrator.setDoSleep(false);

        //instantiate several potentials for selection in combo-box
        P2LennardJones potential = new P2LennardJones(space);
        P2SoftSphericalTruncated p2Truncated = new P2SoftSphericalTruncated(space, potential, 2.5);
        potentialMaster.setPairPotential(species.getLeafType(), species.getLeafType(), p2Truncated);

        //construct box
        Vector dim = space.makeVector();
        dim.E(15);
        box.getBoundary().setBoxSize(dim);
        box.setNMolecules(species, N);
        new ConfigurationLattice(space.D() == 2 ? (new LatticeOrthorhombicHexagonal(space)) : (new LatticeCubicFcc(space)), space).initializeCoordinates(box);
    }

    public IntegratorVelocityVerletFasterer getIntegrator() {
        return integrator;
    }

    public static void main(String[] args) {
        Space space = Space2D.getInstance();
        if (args.length != 0) {
            try {
                int D = Integer.parseInt(args[0]);
                if (D == 3) {
                    space = Space3D.getInstance();
                }
            } catch (NumberFormatException e) {
            }
        }

        LjmdFasterer sim = new LjmdFasterer(space);
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, Long.MAX_VALUE));
    }//end of main

}
