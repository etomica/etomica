/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.mu;


import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.config.ConfigurationLatticeRandom;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorHardFasterer;
import etomica.integrator.IntegratorMDFasterer;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.nbr.list.NeighborListManagerFastererHard;
import etomica.potential.BondingInfo;
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

public class MuFasterer extends Simulation {

    public final SpeciesGeneral speciesA, speciesB;
    public final Box box;
    public final IntegratorHardFasterer integrator;
    public final PotentialComputePairGeneral potentialMasterMu;

    public final P2SquareWellOneSideFasterer potentialAA, potentialAB, potentialBB;
    public final P2HardGeneric potentialAAmu, potentialABmu, potentialBBmu;
    public final P1MagicWallFasterer p1Wall;
    public final ConfigurationLatticeRandom configuration;

    public MuFasterer(Space _space) {
        super(_space);

        //species
        speciesA = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this), true);
        addSpecies(speciesA);
        speciesB = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this), true);
        addSpecies(speciesB);

        box = this.makeBox(new BoundaryRectangularSlit(0, space));

        Vector dim = space.makeVector();
        double xdim = space.D() == 3 ? 30 : 50;
        dim.E(0.5 * xdim);
        dim.setX(0, xdim);
        box.getBoundary().setBoxSize(dim);

        NeighborListManagerFastererHard neighborManager = new NeighborListManagerFastererHard(getSpeciesManager(), box, 2, 4, BondingInfo.noBonding());
        PotentialComputePairGeneral potentialMaster = new PotentialComputePairGeneral(getSpeciesManager(), box, neighborManager);
        potentialMasterMu = new PotentialComputePairGeneral(getSpeciesManager(), box, neighborManager.getCellManager());
        PotentialComputeField pcField = new PotentialComputeField(getSpeciesManager(), box);

        //instantiate several potentials for selection in combo-box
        double sigma = 1.0;
        double lambda = 1.5;
        double epsilon = 1.0;

        potentialAA = new P2SquareWellOneSideFasterer(sigma, lambda, epsilon);
        potentialMaster.setPairPotential(speciesA.getLeafType(), speciesA.getLeafType(), potentialAA);
        potentialAB = new P2SquareWellOneSideFasterer(sigma, lambda, epsilon);
        potentialMaster.setPairPotential(speciesA.getLeafType(), speciesB.getLeafType(), potentialAB);
        potentialBB = new P2SquareWellOneSideFasterer(sigma, lambda, epsilon);
        potentialMaster.setPairPotential(speciesB.getLeafType(), speciesB.getLeafType(), potentialBB);

        potentialAAmu = P2SquareWell.makePotential(sigma, lambda, epsilon);
        potentialMasterMu.setPairPotential(speciesA.getLeafType(), speciesA.getLeafType(), potentialAAmu);
        potentialABmu = P2SquareWell.makePotential(sigma, lambda, epsilon);
        potentialMasterMu.setPairPotential(speciesA.getLeafType(), speciesB.getLeafType(), potentialABmu);
        potentialBBmu = P2SquareWell.makePotential(sigma, lambda, epsilon);
        potentialMasterMu.setPairPotential(speciesB.getLeafType(), speciesB.getLeafType(), potentialBBmu);

        neighborManager.setNeighborRange(4);

        double L = box.getBoundary().getBoxSize().getX(0);
        p1Wall = new P1MagicWallFasterer(L, potentialMaster, neighborManager);
        pcField.setFieldPotential(speciesA.getLeafType(), p1Wall);
        pcField.setFieldPotential(speciesB.getLeafType(), p1Wall);

        int N = 300;  //number of atoms

        //controller and integrator
        integrator = new IntegratorHardFasterer(IntegratorHardFasterer.extractHardPotentials(potentialMaster), IntegratorHardFasterer.extractFieldPotentials(pcField), neighborManager, random, 0.01, 1, box, null, null);
        integrator.setTemperature(1);
        integrator.setIsothermal(true);
        integrator.setThermostat(IntegratorMDFasterer.ThermostatType.ANDERSEN_SINGLE);
        integrator.setThermostatInterval(1);
        getController().addActivity(new ActivityIntegrate(integrator, true));

        //construct box
        box.setNMolecules(speciesA, N);
        box.setNMolecules(speciesB, N);
        configuration = new ConfigurationLatticeRandom(space.D() == 3 ? new LatticeCubicFcc(space) : new LatticeOrthorhombicHexagonal(space), space, random);
        configuration.initializeCoordinates(box);
    }

    public static void main(String[] args) {
        Space space = Space3D.getInstance();
        MuFasterer sim = new MuFasterer(space);
        SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.GRAPHIC_ONLY, "", 1);
        simGraphic.makeAndDisplayFrame();
    }
}
