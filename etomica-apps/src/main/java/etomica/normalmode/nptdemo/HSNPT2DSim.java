/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode.nptdemo;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.integrator.IntegratorHard;
import etomica.lattice.crystal.BasisOrthorhombicHexagonal;
import etomica.lattice.crystal.PrimitiveOrthorhombicHexagonal;
import etomica.nbr.list.NeighborListManagerHard;
import etomica.normalmode.CoordinateDefinition;
import etomica.normalmode.CoordinateDefinitionLeaf;
import etomica.potential.BondingInfo;
import etomica.potential.P2HardGeneric;
import etomica.potential.P2HardSphere;
import etomica.potential.compute.PotentialComputePair;
import etomica.simulation.Simulation;
import etomica.space.Vector;
import etomica.space2d.Space2D;
import etomica.species.SpeciesGeneral;

/**
 * Simple hard-sphere molecular dynamics simulation in 2D.
 *
 * @author David Kofke
 */
 
public class HSNPT2DSim extends Simulation {

    public IntegratorHard integrator;
    public SpeciesGeneral species1;
    public Box box;
    public P2HardGeneric potential;
    public PotentialComputePair potentialMaster;
    public NeighborListManagerHard neighborManager;
    public CoordinateDefinition coordinateDefinition;
    public double pressure;

    public HSNPT2DSim() {
        super(Space2D.getInstance());

        species1 = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this), true);
        addSpecies(species1);
        box = this.makeBox();

        double sigma = 1;
        double neighborRangeFac = 1.6;
        double nbrRange = sigma * neighborRangeFac;
        neighborManager = new NeighborListManagerHard(getSpeciesManager(), box, 2, nbrRange, BondingInfo.noBonding());
        potentialMaster = new PotentialComputePair(getSpeciesManager(), box, neighborManager);

        AtomType leafType1 = species1.getLeafType();
        potential = P2HardSphere.makePotential(sigma);

        potentialMaster.setPairPotential(leafType1, leafType1, potential);

        integrator = new IntegratorHard(potentialMaster.getPairPotentials(), neighborManager, random, 0.05, 1.0, box, getSpeciesManager());
        integrator.setIsothermal(false);
        integrator.setTimeStep(0.05);

        this.getController().addActivity(new ActivityIntegrate(integrator));
        int nx = 10;
        int ny = 6;
        double rho = 1.0;
        box.setNMolecules(species1, nx * ny * 2);
        double bx = 1;
        double by = Math.sqrt(3);
        double v1 = Math.sqrt(3);
        double v2 = 2 / rho;
        bx *= nx * Math.sqrt(v2 / v1);
        by *= ny * Math.sqrt(v2 / v1);
        box.getBoundary().setBoxSize(Vector.of(new double[]{bx, by}));

        coordinateDefinition = new CoordinateDefinitionLeaf(box, new PrimitiveOrthorhombicHexagonal(space, bx / nx), new BasisOrthorhombicHexagonal(), space);
        coordinateDefinition.initializeCoordinates(new int[]{nx, ny});

        potentialMaster.init();
        integrator.getEventManager().removeListener(neighborManager);
    }

}
