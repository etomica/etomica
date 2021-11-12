/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.adsorption;


import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.chem.elements.ElementSimple;
import etomica.config.ConfigurationLattice;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorHard;
import etomica.integrator.IntegratorMC;
import etomica.integrator.IntegratorMD;
import etomica.lattice.LatticeCubicFcc;
import etomica.nbr.cell.NeighborCellManager;
import etomica.nbr.list.NeighborListManagerHard;
import etomica.potential.BondingInfo;
import etomica.potential.P2HardGeneric;
import etomica.potential.P2SquareWell;
import etomica.potential.compute.PotentialComputeField;
import etomica.potential.compute.PotentialComputePair;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularSlit;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesGeneral;

/**
 * Simulation for Adsorption module.
 * Design by Lev Gelb
 *
 * @author Andrew Schultz
 */
public class Adsorption extends Simulation {

    public final SpeciesGeneral speciesA, speciesB;
    public final Box box;
    public final IntegratorHard integratorMD;
    public final IntegratorMC integratorMC;

    public final P2HardGeneric p2AA, p2AB, p2BB;
    public final P1Wall p1WallA, p1WallB;
    public final MyMCMove mcMoveIDA, mcMoveIDB;
    public final ConfigurationLattice config;

    public Adsorption() {
        super(Space3D.getInstance());

        //species
        speciesA = SpeciesGeneral.monatomic(space, AtomType.element(new ElementSimple(this)), true);
        addSpecies(speciesA);
        speciesB = SpeciesGeneral.monatomic(space, AtomType.element(new ElementSimple(this)), true);
        addSpecies(speciesB);

        //construct box
        box = this.makeBox(new BoundaryRectangularSlit(1, 20.0, space));
        NeighborListManagerHard neighborManager = new NeighborListManagerHard(getSpeciesManager(), box, 1, 2, BondingInfo.noBonding());
        neighborManager.setDoDownNeighbors(true);
        PotentialComputePair computePair = new PotentialComputePair(getSpeciesManager(), box, neighborManager);

        PotentialComputeField computeWall = new PotentialComputeField(getSpeciesManager(), box);
        NeighborCellManager neighborManagerMC = new NeighborCellManager(getSpeciesManager(), box, 1, BondingInfo.noBonding());

        double sigma = 1;
        double lambda = 1.5;
        double epsilon = 1.0;
        double epsilonWF = 5.0;

        //potentials
        p2AA = P2SquareWell.makePotential(sigma, lambda, epsilon);
        computePair.setPairPotential(speciesA.getLeafType(), speciesA.getLeafType(), p2AA);
        p2AB = P2SquareWell.makePotential(sigma, lambda, epsilon);
        computePair.setPairPotential(speciesA.getLeafType(), speciesB.getLeafType(), p2AB);
        p2BB = P2SquareWell.makePotential(sigma, lambda, epsilon);
        computePair.setPairPotential(speciesB.getLeafType(), speciesB.getLeafType(), p2BB);
        PotentialComputePair computePairMC = new PotentialComputePair(getSpeciesManager(), box, neighborManagerMC, computePair.getPairPotentials());

        double L = 12 * sigma;
        p1WallA = new P1Wall(L, sigma, sigma / 2, epsilonWF);

        p1WallB = new P1Wall(L, sigma, sigma / 2, epsilonWF);

        computeWall.setFieldPotential(speciesA.getLeafType(), p1WallA);
        computeWall.setFieldPotential(speciesB.getLeafType(), p1WallB);

        integratorMD = new IntegratorHard(IntegratorHard.extractHardPotentials(computePair), IntegratorHard.extractFieldPotentials(computeWall), neighborManager, random, 0.005, 1, box, getSpeciesManager(), null);
        integratorMD.setTimeStep(0.005);
        integratorMD.setIsothermal(true);
        integratorMD.setThermostatInterval(50);

        integratorMD.setThermostat(IntegratorMD.ThermostatType.HYBRID_MC);
        integratorMC = new IntegratorMC(computePairMC, random, 1.0, box);
        integratorMD.setIntegratorMC(integratorMC, 1);

        p1WallA.setThermalize(integratorMD, 0.0, random);
        p1WallB.setThermalize(integratorMD, 0.0, random);

        getController().addActivity(new ActivityIntegrate(integratorMD));

        mcMoveIDA = new MyMCMove(integratorMC, random, space, 0.1, sigma, 1);
        mcMoveIDA.setMu(-12);
        integratorMC.getMoveManager().addMCMove(mcMoveIDA);
        mcMoveIDA.setSpecies(speciesA);
        mcMoveIDA.setBox(box);
        integratorMC.getEventManager().addListener(mcMoveIDA);

        mcMoveIDB = new MyMCMove(integratorMC, random, space, 0.1, sigma, 1);
        mcMoveIDB.setMu(-Double.POSITIVE_INFINITY);
        mcMoveIDB.setSpecies(speciesB);
        mcMoveIDB.setBox(box);
        integratorMC.getEventManager().addListener(mcMoveIDB);

        Vector dim = space.makeVector();
        dim.E(8 * sigma);
        dim.setX(1, L);
        box.getBoundary().setBoxSize(dim);

        box.setNMolecules(speciesA, 40);

        config = new ConfigurationLattice(new LatticeCubicFcc(space), space);
        config.setBoundaryPadding(1.2345);
        config.initializeCoordinates(box);
    }

    public static void main(String[] args) {
        Space space = Space3D.getInstance();

        Adsorption sim = new Adsorption();
        SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.GRAPHIC_ONLY, "Catalysis", 1);
        simGraphic.makeAndDisplayFrame();
    }
}
