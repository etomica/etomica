/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.mu;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorHard;
import etomica.integrator.IntegratorMD.ThermostatType;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.modules.chainequilibrium.ConfigurationLatticeRandom;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.P1HardBoundary;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularSlit;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;

public class Mu extends Simulation {
    
    private static final long serialVersionUID = 1L;
    public final SpeciesSpheresMono speciesA, speciesB;
    public final Box box;
    public final IntegratorHard integrator;
    public final ActivityIntegrate activityIntegrate;
    public final P2SquareWellOneSide potentialAA, potentialAB, potentialBB;
    public final P1MagicWall p1WallA, p1WallB;
    public final P1HardBoundary p1BoundaryA, p1BoundaryB;
    public final ConfigurationLatticeRandom configuration;
    
    public Mu(Space _space) {
        super(_space);

        //species
        speciesA = new SpeciesSpheresMono(this, space);
        speciesA.setIsDynamic(true);
        addSpecies(speciesA);
        speciesB = new SpeciesSpheresMono(this, space);
        speciesB.setIsDynamic(true);
        addSpecies(speciesB);

        box = this.makeBox(new BoundaryRectangularSlit(0, space));
        PotentialMasterList potentialMaster = new PotentialMasterList(this, 4, space); //List(this, 2.0);

        int N = 300;  //number of atoms

        double sigma = 1.0;
        double lambda = 1.5;
        double epsilon = 1.0;

        //controller and integrator
        integrator = new IntegratorHard(this, potentialMaster, box);
        integrator.setTimeStep(0.02);
        integrator.setTemperature(1);
        integrator.setIsothermal(true);
        integrator.setThermostat(ThermostatType.ANDERSEN_SINGLE);
        integrator.setThermostatInterval(1);
        activityIntegrate = new ActivityIntegrate(integrator, 0, true);
        getController().addAction(activityIntegrate);

        //instantiate several potentials for selection in combo-box
        potentialAA = new P2SquareWellOneSide(space, sigma, lambda, epsilon, true);
        potentialMaster.addPotential(potentialAA, new AtomType[]{speciesA.getLeafType(), speciesA.getLeafType()});
        potentialAB = new P2SquareWellOneSide(space, sigma, lambda, epsilon, true);
        potentialMaster.addPotential(potentialAB, new AtomType[]{speciesA.getLeafType(), speciesB.getLeafType()});
        potentialBB = new P2SquareWellOneSide(space, sigma, lambda, epsilon, true);
        potentialMaster.addPotential(potentialBB, new AtomType[]{speciesB.getLeafType(), speciesB.getLeafType()});

        p1BoundaryA = new P1HardBoundary(space);
        p1BoundaryA.setActive(0, false, true);
        p1BoundaryA.setActive(0, true, true);
        p1BoundaryA.setActive(1, false, false);
        p1BoundaryA.setActive(1, true, false);
        if (space.D() == 3) {
            p1BoundaryA.setActive(2, false, false);
            p1BoundaryA.setActive(2, true, false);
        }
        p1BoundaryB = new P1HardBoundary(space);
        p1BoundaryB.setActive(0, false, true);
        p1BoundaryB.setActive(0, true, true);
        p1BoundaryB.setActive(1, false, false);
        p1BoundaryB.setActive(1, true, false);
        if (space.D() == 3) {
            p1BoundaryB.setActive(2, false, false);
            p1BoundaryB.setActive(2, true, false);
        }
        potentialMaster.addPotential(p1BoundaryA, new AtomType[]{speciesA.getLeafType()});
        potentialMaster.addPotential(p1BoundaryB, new AtomType[]{speciesB.getLeafType()});
        p1WallA = new P1MagicWall(space, potentialMaster);
        potentialMaster.addPotential(p1WallA, new AtomType[]{speciesA.getLeafType()});
        p1WallB = new P1MagicWall(space, potentialMaster);
        potentialMaster.addPotential(p1WallB, new AtomType[]{speciesB.getLeafType()});
        //construct box
        Vector dim = space.makeVector();
        double xdim = space.D() == 3 ? 30 : 50;
        dim.E(0.5 * xdim);
        dim.setX(0, xdim);
        box.getBoundary().setBoxSize(dim);
        box.setNMolecules(speciesA, N);
        box.setNMolecules(speciesB, N);
        configuration = new ConfigurationLatticeRandom(space.D() == 3 ? new LatticeCubicFcc(space) : new LatticeOrthorhombicHexagonal(space), space, random);
        configuration.initializeCoordinates(box);
        integrator.getEventManager().addListener(potentialMaster.getNeighborManager(box));
    }
    
    public static void main(String[] args) {
        Space space = Space3D.getInstance();
        Mu sim = new Mu(space);
        SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.GRAPHIC_ONLY, "", 1);
        simGraphic.makeAndDisplayFrame();
    }
}
