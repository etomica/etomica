/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.catalysis;


import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.chem.elements.Carbon;
import etomica.chem.elements.ElementSimple;
import etomica.chem.elements.Oxygen;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorHard;
import etomica.integrator.IntegratorMD.ThermostatType;
import etomica.nbr.list.NeighborListManagerHard;
import etomica.potential.BondingInfo;
import etomica.potential.P1HardBoundary;
import etomica.potential.P2HardGeneric;
import etomica.potential.P2SquareWell;
import etomica.potential.compute.PotentialComputeField;
import etomica.potential.compute.PotentialComputePairGeneral;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularSlit;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesGeneral;
import etomica.units.Calorie;
import etomica.units.Kelvin;
import etomica.units.Mole;
import etomica.util.random.IRandom;
import etomica.util.random.RandomMersenneTwister;
import etomica.util.random.RandomNumberGeneratorUnix;

/**
 * Simulation for Catalysis module.
 * Design by Ken Benjamin
 *
 * @author Andrew Schultz
 */
public class Catalysis extends Simulation {

    public final SpeciesGeneral speciesO, speciesC, speciesSurface;
    public final Box box;
    public final IntegratorHard integrator;

    public final P2SquareWellBonding potentialOO;
    public final P2SquareWellBondingCO potentialCO;
    public final P2HardGeneric potentialCC;
    public final P2SquareWellSurface potentialCS, potentialOS;
    public final ConfigurationCatalysis config;
    public final InteractionTracker interactionTracker;
    public final ReactionManagerCO reactionManagerCO;

    public Catalysis(Space space, int nCellsZ) {
        this(space, nCellsZ, new RandomMersenneTwister(RandomNumberGeneratorUnix.getRandSeedArray()));
    }

    public Catalysis(Space _space, int nCellsZ, IRandom random) {
        super(_space);
        this.setRandom(random);
        //species
        speciesO = SpeciesGeneral.monatomic(space, AtomType.element(Oxygen.INSTANCE), true);
        addSpecies(speciesO);
        speciesC = SpeciesGeneral.monatomic(space, AtomType.element(Carbon.INSTANCE), true);
        addSpecies(speciesC);
        speciesSurface = SpeciesGeneral.monatomic(space, AtomType.element(new ElementSimple("Surface", Double.POSITIVE_INFINITY)), true);
        addSpecies(speciesSurface);

        //construct box
        box = this.makeBox(new BoundaryRectangularSlit(1, 20.0, space));

        NeighborListManagerHard neighborManager = new NeighborListManagerHard(getSpeciesManager(), box, 2, 9, BondingInfo.noBonding());
        neighborManager.setDoDownNeighbors(true);
        PotentialComputePairGeneral potentialMaster = new PotentialComputePairGeneral(getSpeciesManager(), box, neighborManager);
        //controller and integrator

        double sigmaO = 3.6;
        double sigmaC = 3.8;
        double sigmaS = 3.7;
        double epsilonO = 1000 * Mole.UNIT.fromSim(Calorie.UNIT.toSim(0.15));
        double epsilonC = 1000 * Mole.UNIT.fromSim(Calorie.UNIT.toSim(0.08));
        double epsilonOS = Kelvin.UNIT.toSim(500);
        double epsilonCS = Kelvin.UNIT.toSim(2000);

        //potentials
        interactionTracker = new InteractionTracker(box, speciesSurface);

        int minOSites = 2, minCSites = 2;

        potentialOO = new P2SquareWellBonding(interactionTracker.getAgentManager(), sigmaO, 1.3, epsilonO, minOSites, Kelvin.UNIT.toSim(200), Kelvin.UNIT.toSim(400), 7.4);
        potentialMaster.setPairPotential(speciesO.getLeafType(), speciesO.getLeafType(), potentialOO);

        potentialCO = new P2SquareWellBondingCO(space, interactionTracker.getAgentManager(), 0.5 * (sigmaO + sigmaC), 1.1, Math.sqrt(epsilonC * epsilonO), 20, Kelvin.UNIT.toSim(400), Kelvin.UNIT.toSim(7500), 7.4, box.getBoundary());
        potentialMaster.setPairPotential(speciesO.getLeafType(), speciesC.getLeafType(), potentialCO);

        potentialCC = P2SquareWell.makePotential(sigmaC, 1.3, epsilonC);
        potentialMaster.setPairPotential(speciesC.getLeafType(), speciesC.getLeafType(), potentialCC);

        potentialOS = new P2SquareWellSurface(interactionTracker.getAgentManager(), 0.5 * (sigmaO + sigmaS), 1.3, epsilonOS, minOSites);
        potentialMaster.setPairPotential(speciesO.getLeafType(), speciesSurface.getLeafType(), potentialOS);

        potentialCS = new P2SquareWellSurface(interactionTracker.getAgentManager(), 0.5 * (sigmaC + sigmaS), 1.3, epsilonCS, minCSites);
        potentialMaster.setPairPotential(speciesC.getLeafType(), speciesSurface.getLeafType(), potentialCS);

        potentialCO.setCSPotential(potentialCS);

        PotentialComputeField pcField = new PotentialComputeField(getSpeciesManager(), box);
        P1HardBoundary p1HardWallO = new P1HardBoundary(space, true, box);
        p1HardWallO.setActive(0, true, false);
        p1HardWallO.setActive(0, false, false);
        p1HardWallO.setActive(1, true, false);
        p1HardWallO.setActive(1, false, true);
        p1HardWallO.setActive(2, true, false);
        p1HardWallO.setActive(2, false, false);
        pcField.setFieldPotential(speciesO.getLeafType(), p1HardWallO);
        P1HardBoundary p1HardWallC = new P1HardBoundary(space, true, box);
        p1HardWallC.setActive(0, true, false);
        p1HardWallC.setActive(0, false, false);
        p1HardWallC.setActive(1, true, false);
        p1HardWallC.setActive(1, false, true);
        p1HardWallC.setActive(2, true, false);
        p1HardWallC.setActive(2, false, false);
        pcField.setFieldPotential(speciesC.getLeafType(), p1HardWallC);

        integrator = new IntegratorHard(IntegratorHard.extractHardPotentials(potentialMaster), IntegratorHard.extractFieldPotentials(pcField),
                neighborManager, random, 0.005, Kelvin.UNIT.toSim(600), box, getSpeciesManager(), null);
        integrator.setIsothermal(true);
        integrator.setThermostat(ThermostatType.ANDERSEN_SINGLE);
        integrator.setThermostatInterval(1);
        getController().addActivity(new ActivityIntegrate(integrator));
        reactionManagerCO = new ReactionManagerCO(this);
        reactionManagerCO.setnReactCO(minCSites);
        integrator.getEventManager().addListener(reactionManagerCO);

        Vector dim = space.makeVector();
        dim.E(100);
        dim.setX(1, 60);
        box.getBoundary().setBoxSize(dim);
        int nCO = 40, nO2 = 40;

        config = new ConfigurationCatalysis(this, space, speciesSurface, speciesC, speciesO, interactionTracker.getAgentManager());
        config.setNCellsX(nCellsZ * 3 / 2);
        config.setNCellsZ(nCellsZ);
        config.setCellSizeX(sigmaS);
        config.setCellSizeZ(sigmaS * Math.sqrt(3));
        config.setNumCO(nCO);
        config.setNumO2(nO2);
        config.initializeCoordinates(box);
    }

    public static void main(String[] args) {
        Space space = Space3D.getInstance();
        int nCellsZ = 20;

        Catalysis sim = new Catalysis(space, nCellsZ);
        SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.GRAPHIC_ONLY, "Catalysis", 1);
        simGraphic.makeAndDisplayFrame();
    }
}
