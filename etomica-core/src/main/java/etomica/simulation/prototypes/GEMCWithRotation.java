/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.simulation.prototypes;

import etomica.action.BoxImposePbc;
import etomica.action.BoxInflate;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.integrator.IntegratorGEMC;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.integrator.mcmove.MCMoveManager;
import etomica.integrator.mcmove.MCMoveRotate;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.lattice.SpaceLattice;
import etomica.listener.IntegratorListenerAction;
import etomica.potential.P2LennardJones;
import etomica.potential.PotentialMaster;
import etomica.potential.PotentialMasterMonatomic;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.species.SpeciesSpheresRotating;

/**
 * Simple Gibbs-ensemble Monte Carlo simulation of rotating molecules.
 */
//in present form uses just a LJ potential, so orientation is irrelevant
public class GEMCWithRotation extends Simulation {

    private static final long serialVersionUID = 1L;
    public Box box1, box2;
    public IntegratorGEMC integrator;
    public SpeciesSpheresRotating species;
    public P2LennardJones potential;

    public GEMCWithRotation(Space _space) {
        super(_space);
        double sigma = 1.2;
        PotentialMaster potentialMaster = new PotentialMasterMonatomic(this);
        integrator = new IntegratorGEMC(getRandom(), space);
        integrator.setEventInterval(400);
        ActivityIntegrate activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);
        activityIntegrate.setSleepPeriod(1);

        species = new SpeciesSpheresRotating(this, space);
        addSpecies(species);

        box1 = new Box(space);
        addBox(box1);
        box1.setNMolecules(species, 200);

        IntegratorMC integratorMC1 = new IntegratorMC(this, potentialMaster);
        integratorMC1.setBox(box1);
        integratorMC1.setTemperature(0.420);
        MCMoveManager moveManager = integratorMC1.getMoveManager();
        moveManager.addMCMove(new MCMoveRotate(potentialMaster, getRandom(), space));
        moveManager.addMCMove(new MCMoveAtom(random, potentialMaster, space));
        integrator.addIntegrator(integratorMC1);


        box2 = new Box(space);
        addBox(box2);
        box2.setNMolecules(species, 200);
        IntegratorMC integratorMC2 = new IntegratorMC(this, potentialMaster);
        integratorMC2.setBox(box2);
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

        potential = new P2LennardJones(space);
        potential.setSigma(sigma);

        potentialMaster.addPotential(potential, new AtomType[]{species.getLeafType(), species.getLeafType()});

        integratorMC1.getEventManager().addListener(new IntegratorListenerAction(new BoxImposePbc(box1, space)));
        integratorMC2.getEventManager().addListener(new IntegratorListenerAction(new BoxImposePbc(box2, space)));

        BoxInflate inflater = new BoxInflate(box2, space);
        inflater.setTargetDensity(0.1);
        inflater.actionPerformed();
    }
}
