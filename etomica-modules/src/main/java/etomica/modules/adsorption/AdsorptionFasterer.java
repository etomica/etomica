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
import etomica.integrator.IntegratorHardFasterer;
import etomica.integrator.IntegratorMCFasterer;
import etomica.integrator.IntegratorMDFasterer;
import etomica.lattice.LatticeCubicFcc;
import etomica.nbr.cell.NeighborCellManagerFasterer;
import etomica.nbr.list.NeighborListManagerFastererHard;
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
import etomica.util.random.RandomMersenneTwister;

/**
 * Simulation for Adsorption module.
 * Design by Lev Gelb
 *
 * @author Andrew Schultz
 */
public class AdsorptionFasterer extends Simulation {

    public final SpeciesGeneral speciesA, speciesB;
    public final Box box;
    public final IntegratorHardFasterer integratorMD;
    public final IntegratorMCFasterer integratorMC;

    public final P2HardGeneric p2AA, p2AB, p2BB;
    public final P1WallFasterer p1WallA, p1WallB;
    public final MyMCMoveFasterer mcMoveIDA, mcMoveIDB;
    public final ConfigurationLattice config;

    public AdsorptionFasterer() {
        super(Space3D.getInstance());
        setRandom(new RandomMersenneTwister(1));

        //species
        speciesA = SpeciesGeneral.monatomic(space, AtomType.element(new ElementSimple(this)), true);
        addSpecies(speciesA);
        speciesB = SpeciesGeneral.monatomic(space, AtomType.element(new ElementSimple(this)), true);
        addSpecies(speciesB);

        //construct box
        box = this.makeBox(new BoundaryRectangularSlit(1, 20.0, space));
        NeighborListManagerFastererHard neighborManager = new NeighborListManagerFastererHard(getSpeciesManager(), box, 1, 2, BondingInfo.noBonding());
        neighborManager.setDoDownNeighbors(true);
        PotentialComputePair computePair = new PotentialComputePair(this, box, neighborManager);

        PotentialComputeField computeWall = new PotentialComputeField(this, box);

        //controller and integrator
        integratorMD = new IntegratorHardFasterer(computePair, computeWall, neighborManager, random, 0.005, 1, box, null);
        integratorMD.setTimeStep(0.005);
        integratorMD.setIsothermal(true);
        integratorMD.setThermostatInterval(50);

        integratorMD.setThermostat(IntegratorMDFasterer.ThermostatType.HYBRID_MC);
        NeighborCellManagerFasterer neighborManagerMC = new NeighborCellManagerFasterer(getSpeciesManager(), box, 1, BondingInfo.noBonding());
        PotentialComputePair computePairMC = new PotentialComputePair(this, box, neighborManagerMC);
        integratorMC = new IntegratorMCFasterer(computePairMC, random, 1.0, box);
        integratorMD.setIntegratorMC(integratorMC, 1);

        getController().addActivity(new ActivityIntegrate(integratorMD));

        double sigma = 1;
        double lambda = 1.5;
        double epsilon = 1.0;
        double epsilonWF = 5.0;

        //potentials
        mcMoveIDA = new MyMCMoveFasterer(integratorMC, random, space, 0.1, sigma, 1);
        mcMoveIDA.setMu(-12);
        integratorMC.getMoveManager().addMCMove(mcMoveIDA);
        mcMoveIDA.setSpecies(speciesA);
        mcMoveIDA.setBox(box);
        integratorMC.getEventManager().addListener(mcMoveIDA);

        mcMoveIDB = new MyMCMoveFasterer(integratorMC, random, space, 0.1, sigma, 1);
        mcMoveIDB.setMu(-Double.POSITIVE_INFINITY);
        mcMoveIDB.setSpecies(speciesB);
        mcMoveIDB.setBox(box);
        integratorMC.getEventManager().addListener(mcMoveIDB);


        p2AA = P2SquareWell.makePotential(sigma, lambda, epsilon);
        computePair.setPairPotential(speciesA.getLeafType(), speciesA.getLeafType(), p2AA);
        p2AB = P2SquareWell.makePotential(sigma, lambda, epsilon);
        computePair.setPairPotential(speciesA.getLeafType(), speciesB.getLeafType(), p2AB);
        p2BB = P2SquareWell.makePotential(sigma, lambda, epsilon);
        computePair.setPairPotential(speciesB.getLeafType(), speciesB.getLeafType(), p2BB);
        computePairMC.setPairPotentials(computePair.getPairPotentials());

        double L = 12 * sigma;
        p1WallA = new P1WallFasterer(L, sigma, sigma / 2, epsilonWF);
        p1WallA.setThermalize(integratorMD, 0.0, random);

        p1WallB = new P1WallFasterer(L, sigma, sigma / 2, epsilonWF);
        p1WallB.setThermalize(integratorMD, 0.0, random);

        computeWall.setFieldPotential(speciesA.getLeafType(), p1WallA);
        computeWall.setFieldPotential(speciesB.getLeafType(), p1WallB);

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

        AdsorptionFasterer sim = new AdsorptionFasterer();
        SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.GRAPHIC_ONLY, "Catalysis", 1);
        simGraphic.makeAndDisplayFrame();
    }
}
