/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.simulation.prototypes;

import etomica.action.BoxImposePbc;
import etomica.action.BoxInflate;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.AtomTypeOriented;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.graphics.DisplayBoxCanvasG3DSys;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorGEMC;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.integrator.mcmove.MCMoveManager;
import etomica.integrator.mcmove.MCMoveRotate;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.lattice.SpaceLattice;
import etomica.potential.P2HardAssociationCone;
import etomica.potential.PotentialMaster;
import etomica.potential.PotentialMasterMonatomic;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresRotating;

import java.awt.*;

/**
 * Simple Gibbs-ensemble Monte Carlo simulation of rotating molecules.
 * <p>
 * Molecules interact with a hard sphere potential and a directional square well association.
 */
public class GEMCWithRotation extends Simulation {

    public Box box1, box2;
    public IntegratorGEMC integrator;
    public SpeciesSpheresRotating species;
    public P2HardAssociationCone potential;

    public GEMCWithRotation() {
        this(Space3D.getInstance());
    }

    public GEMCWithRotation(Space _space) {
        super(_space);

        species = new SpeciesSpheresRotating(this, space);
        addSpecies(species);

        double sigma = 1.2;
        PotentialMaster potentialMaster = new PotentialMasterMonatomic(this);
        integrator = new IntegratorGEMC(getRandom(), space);
        integrator.setTemperature(0.420);
        integrator.setEventInterval(400);

        getController().addActivity(new ActivityIntegrate(integrator)).setSleepPeriod(1);

        box1 = this.makeBox();
        box1.setNMolecules(species, 200);

        IntegratorMC integratorMC1 = new IntegratorMC(this, potentialMaster, box1);
        integratorMC1.setTemperature(0.420);
        MCMoveManager moveManager = integratorMC1.getMoveManager();
        moveManager.addMCMove(new MCMoveRotate(potentialMaster, getRandom(), space));
        moveManager.addMCMove(new MCMoveAtom(random, potentialMaster, space));
        integrator.addIntegrator(integratorMC1);


        box2 = this.makeBox();
        box2.setNMolecules(species, 200);
        IntegratorMC integratorMC2 = new IntegratorMC(this, potentialMaster, box2);
        integratorMC2.setTemperature(0.420);
        moveManager = integratorMC2.getMoveManager();
        moveManager.addMCMove(new MCMoveRotate(potentialMaster, getRandom(), space));
        moveManager.addMCMove(new MCMoveAtom(random, potentialMaster, space));
        // GEMC integrator adds volume and molecule exchange moves once
        // it has 2 integrators
        integrator.addIntegrator(integratorMC2);

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

        potentialMaster.addPotential(potential, new AtomType[]{species.getLeafType(), species.getLeafType()});

        integratorMC1.getEventManager().addListener(new IntegratorListenerAction(new BoxImposePbc(box1, space)));
        integratorMC2.getEventManager().addListener(new IntegratorListenerAction(new BoxImposePbc(box2, space)));

        BoxInflate inflater = new BoxInflate(box2, space);
        inflater.setTargetDensity(0.1);
        inflater.actionPerformed();
    }

    public static void main(String[] args) {

        double sigma = 1.2;
        final GEMCWithRotation sim = new GEMCWithRotation();
        final SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
        double conedia = sigma * sim.potential.getTheta();
        DisplayBoxCanvasG3DSys.OrientedSite[] os = new DisplayBoxCanvasG3DSys.OrientedSite[]{new DisplayBoxCanvasG3DSys.OrientedSite(sigma / 2, Color.BLUE, conedia)};
        ((DisplayBoxCanvasG3DSys) (simGraphic.getDisplayBox(sim.box1).canvas)).setOrientationSites(((AtomTypeOriented) sim.species.getLeafType()), os);
        ((DisplayBoxCanvasG3DSys) (simGraphic.getDisplayBox(sim.box2).canvas)).setOrientationSites(((AtomTypeOriented) sim.species.getLeafType()), os);
        simGraphic.makeAndDisplayFrame();
    }
}
