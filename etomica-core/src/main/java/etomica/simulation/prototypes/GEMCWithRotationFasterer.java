/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.simulation.prototypes;

import etomica.action.BoxImposePbc;
import etomica.action.BoxInflate;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomTypeOriented;
import etomica.box.Box;
import etomica.chem.elements.ElementSimple;
import etomica.config.ConfigurationLattice;
import etomica.graphics.DisplayBoxCanvasG3DSys;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorGEMC;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorMCFasterer;
import etomica.integrator.IntegratorManagerMC;
import etomica.integrator.mcmove.MCMoveAtomFasterer;
import etomica.integrator.mcmove.MCMoveAtomRotateFasterer;
import etomica.integrator.mcmove.MCMoveManager;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.lattice.SpaceLattice;
import etomica.potential.P2HardAssociationCone;
import etomica.potential.compute.NeighborManagerSimple;
import etomica.potential.compute.PotentialComputePairGeneral;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.SpeciesGeneral;
import etomica.species.SpeciesSpheresRotating;

import java.awt.*;

/**
 * Simple Gibbs-ensemble Monte Carlo simulation of rotating molecules.
 * <p>
 * Molecules interact with a hard sphere potential and a directional square well association.
 */
public class GEMCWithRotationFasterer extends Simulation {

    public Box box1, box2;
    public IntegratorManagerMC integrator;
    public SpeciesGeneral species;
    public P2HardAssociationCone potential;

    public GEMCWithRotationFasterer() {
        this(Space3D.getInstance());
    }

    public GEMCWithRotationFasterer(Space _space) {
        super(_space);

        species = SpeciesSpheresRotating.create(space, new ElementSimple(this), false, true);
        addSpecies(species);

        box1 = this.makeBox();
        box2 = this.makeBox();

        double sigma = 1.2;
        NeighborManagerSimple neighborManager1 = new NeighborManagerSimple(box1);
        NeighborManagerSimple neighborManager2 = new NeighborManagerSimple(box2);
        PotentialComputePairGeneral potentialMaster1 = new PotentialComputePairGeneral(getSpeciesManager(), box1, neighborManager1);
        PotentialComputePairGeneral potentialMaster2 = new PotentialComputePairGeneral(getSpeciesManager(), box2, neighborManager2);

//        getController().setSleepPeriod(1);

        box1.setNMolecules(species, 200);

        IntegratorMCFasterer integratorMC1 = new IntegratorMCFasterer(potentialMaster1, this.getRandom(), 1.0, box1);
        integratorMC1.setTemperature(0.420);
        MCMoveManager moveManager = integratorMC1.getMoveManager();
        moveManager.addMCMove(new MCMoveAtomRotateFasterer(random, potentialMaster1, box1));
        moveManager.addMCMove(new MCMoveAtomFasterer(random, potentialMaster1, box1));

        box2.setNMolecules(species, 200);
        IntegratorMCFasterer integratorMC2 = new IntegratorMCFasterer(potentialMaster2, this.getRandom(), 1.0, box2);
        integratorMC2.setTemperature(0.420);
        moveManager = integratorMC2.getMoveManager();
        moveManager.addMCMove(new MCMoveAtomRotateFasterer(random, potentialMaster2, box2));
        moveManager.addMCMove(new MCMoveAtomFasterer(random, potentialMaster2, box2));

        // GEMC integrator adds volume and molecule exchange moves
        integrator = IntegratorGEMC.buildGEMC(integratorMC1, integratorMC2, getRandom(), space);
        integrator.setTemperature(0.420);
        integrator.setEventInterval(400);
        getController().addActivity(new ActivityIntegrate(integrator));

        SpaceLattice lattice;
        if (space.D() == 2) {
            lattice = new LatticeOrthorhombicHexagonal(space);
        } else {
            lattice = new LatticeCubicFcc(space);
        }
        ConfigurationLattice config = new ConfigurationLattice(lattice, space);
        config.initializeCoordinates(box1);
        config.initializeCoordinates(box2);

        potential = new P2HardAssociationCone(space);
        potential.setSigma(sigma);
        potential.setWellCutoffFactor(1.5);

        potentialMaster1.setPairPotential(species.getLeafType(), species.getLeafType(), potential);
        potentialMaster2.setPairPotential(species.getLeafType(), species.getLeafType(), potential);

        integratorMC1.getEventManager().addListener(new IntegratorListenerAction(new BoxImposePbc(box1, space)));
        integratorMC2.getEventManager().addListener(new IntegratorListenerAction(new BoxImposePbc(box2, space)));

        BoxInflate inflater = new BoxInflate(box2, space);
        inflater.setTargetDensity(0.1);
        inflater.actionPerformed();
    }

    public static void main(String[] args) {

        double sigma = 1.2;
        final GEMCWithRotationFasterer sim = new GEMCWithRotationFasterer();
        final SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
        double conedia = sigma * sim.potential.getTheta();
        DisplayBoxCanvasG3DSys.OrientedSite[] os = new DisplayBoxCanvasG3DSys.OrientedSite[]{new DisplayBoxCanvasG3DSys.OrientedSite(sigma / 2, Color.BLUE, conedia)};
        ((DisplayBoxCanvasG3DSys) (simGraphic.getDisplayBox(sim.box1).canvas)).setOrientationSites(((AtomTypeOriented) sim.species.getLeafType()), os);
        ((DisplayBoxCanvasG3DSys) (simGraphic.getDisplayBox(sim.box2).canvas)).setOrientationSites(((AtomTypeOriented) sim.species.getLeafType()), os);
        simGraphic.makeAndDisplayFrame();
    }
}
