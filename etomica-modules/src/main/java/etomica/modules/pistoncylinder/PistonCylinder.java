/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.pistoncylinder;


import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.chem.elements.ElementSimple;
import etomica.config.ConfigurationLattice;
import etomica.integrator.IntegratorMD;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.potential.P2HardGeneric;
import etomica.potential.compute.NeighborManagerSimpleHard;
import etomica.potential.compute.PotentialComputeField;
import etomica.potential.compute.PotentialComputePairGeneral;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space2d.Vector2D;
import etomica.space3d.Vector3D;
import etomica.species.SpeciesGeneral;
import etomica.units.Bar;

/**
 * Simple hard-sphere MD in piston-cylinder apparatus
 */
public class PistonCylinder extends Simulation {

    private final int INIT_NUM_MOLECULES = 100;
    public IntegratorHardPiston integrator;
    public P1HardMovingBoundary pistonPotential;
    public SpeciesGeneral species;
    public Box box;

    public double lambda;
    public ConfigurationLattice config;

    public PistonCylinder(int D) {
        super(Space.getInstance(D));
        species = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this), true);
        ((ElementSimple) species.getLeafType().getElement()).setMass(16);
        addSpecies(species);

        lambda = 2.0;
        double sigma = 4.0;
        box = this.makeBox(new BoundaryPistonCylinder(space));
        box.setNMolecules(species, INIT_NUM_MOLECULES);
        Vector newDim;
        if (space.D() == 2) {
            config = new ConfigurationLattice(new LatticeOrthorhombicHexagonal(space), space);
            newDim = new Vector2D(80, 150);
        } else {
            config = new ConfigurationLattice(new LatticeCubicFcc(space), space);
            newDim = new Vector3D(30, 80, 30);
        }
        config.setBoundaryPadding(sigma);
        box.getBoundary().setBoxSize(newDim);
        config.initializeCoordinates(box);

        NeighborManagerSimpleHard neighborManager = new NeighborManagerSimpleHard(box);
        PotentialComputePairGeneral potentialMaster = new PotentialComputePairGeneral(getSpeciesManager(), box, neighborManager);

        P2HardGeneric potentialSW = new P2HardGeneric(new double[]{sigma, sigma * lambda}, new double[]{Double.POSITIVE_INFINITY, -31.875}, true);
        potentialMaster.setPairPotential(species.getLeafType(), species.getLeafType(), potentialSW);

        PotentialComputeField pcField = new PotentialComputeField(getSpeciesManager(), box);

        pistonPotential = new P1HardMovingBoundary(box, 1, 400);
        pistonPotential.setCollisionRadius(sigma * 0.5);
        pistonPotential.setActive(0, true, true);  // left wall
        pistonPotential.setActive(0, false, true); // right wall
        if (D == 3) {
            pistonPotential.setActive(1, true, true); // bottom wall
            pistonPotential.setActive(1, false, false); // top wall
            pistonPotential.setActive(2, true, true);  // front wall
            pistonPotential.setActive(2, false, true); // back wall

            pistonPotential.setWallPosition(box.getBoundary().getBoxSize().getX(1) * 0.5);
            pistonPotential.setWallVelocity(-0.5);
            pistonPotential.setPressure(-Bar.UNIT.toSim(1.0));
        } else {
            pistonPotential.setActive(1, true, false); // top wall
            pistonPotential.setActive(1, false, true); // bottom wall

            pistonPotential.setWallPosition(-box.getBoundary().getBoxSize().getX(1) * 0.5);
            pistonPotential.setWallVelocity(0.5);
            pistonPotential.setPressure(Bar.UNIT.toSim(100.0));
        }
        pistonPotential.setThickness(1.0);
        pcField.setFieldPotential(species.getLeafType(), pistonPotential);
//        potentialMaster.addPotential(pistonPotential, new AtomType[]{species.getLeafType()});
        ((BoundaryPistonCylinder) box.getBoundary()).setPistonPotential(pistonPotential);

        integrator = new IntegratorHardPiston(potentialMaster.getPairPotentials(), pcField.getFieldPotentials(),
                neighborManager, random, 1.0, 1.0, box, pistonPotential, getSpeciesManager());
        integrator.setIsothermal(true);
        integrator.setThermostatInterval(1);
        integrator.setThermostat(IntegratorMD.ThermostatType.ANDERSEN_SINGLE);
        integrator.setTimeStep(1.0);
        getController().addActivity(new ActivityIntegrate(integrator, true));

    }
}
