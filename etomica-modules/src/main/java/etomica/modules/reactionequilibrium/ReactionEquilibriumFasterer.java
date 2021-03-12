/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.reactionequilibrium;

import etomica.action.BoxImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomLeafAgentManager.AgentSource;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.config.Configuration;
import etomica.config.ConfigurationLattice;
import etomica.data.meter.MeterTemperature;
import etomica.integrator.IntegratorHardFasterer;
import etomica.integrator.IntegratorListenerAction;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.potential.IPotentialPair;
import etomica.potential.compute.NeighborIterator;
import etomica.potential.compute.NeighborManagerSimpleHard;
import etomica.potential.compute.PotentialComputePairGeneral;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Vector;
import etomica.space2d.Space2D;
import etomica.species.SpeciesGeneral;

public class ReactionEquilibriumFasterer extends Simulation implements AgentSource<IAtom> {

    public IntegratorHardFasterer integratorHard1;
    public java.awt.Component display;
    public Box box;
    public MeterTemperature thermometer;
    public SpeciesGeneral speciesA;
    public SpeciesGeneral speciesB;
    public P2SquareWellBondedFasterer AAbonded;
    public P2SquareWellBondedFasterer ABbonded;
    public P2SquareWellBondedFasterer BBbonded;
    public MeterDimerFraction meterDimerFraction;
    public PotentialComputePairGeneral potentialMaster;
    public NeighborManagerSimpleHard neighborManager;

    private final AtomLeafAgentManager<IAtom> agentManager;

    public ReactionEquilibriumFasterer() {
        super(Space2D.getInstance());

        speciesA = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this), true);
        speciesB = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this), true);
        addSpecies(speciesA);
        addSpecies(speciesB);

        box = this.makeBox(new BoundaryRectangularPeriodic(space, 30.0));
        box.setNMolecules(speciesA, 30);
        box.setNMolecules(speciesB, 30);

        neighborManager = new NeighborManagerSimpleHard(box);
        potentialMaster = new PotentialComputePairGeneral(getSpeciesManager(), box, neighborManager);

        double diameter = 1.0;

        agentManager = new AtomLeafAgentManager<>(this, box);

        //potentials
        AAbonded = new P2SquareWellBondedFasterer(agentManager, 0.5 * diameter, //core
                2.0, //well multiplier
                1.0);
        ABbonded = new P2SquareWellBondedFasterer(agentManager, 0.5 * diameter, //core
                2.0, //well multiplier
                1.0);
        BBbonded = new P2SquareWellBondedFasterer(agentManager, 0.5 * diameter, //core
                2.0, //well multiplier
                1.0);
        potentialMaster.setPairPotential(speciesA.getLeafType(), speciesA.getLeafType(), AAbonded);
        potentialMaster.setPairPotential(speciesA.getLeafType(), speciesB.getLeafType(), ABbonded);
        potentialMaster.setPairPotential(speciesB.getLeafType(), speciesB.getLeafType(), BBbonded);

        integratorHard1 = new IntegratorHardFasterer(IntegratorHardFasterer.extractHardPotentials(potentialMaster),
                neighborManager, random, 0.05, 1.0, box);
        integratorHard1.setMaxCollisionDiameter(speciesA.getLeafType(), diameter);
        integratorHard1.setMaxCollisionDiameter(speciesB.getLeafType(), diameter);
        integratorHard1.setIsothermal(true);

        meterDimerFraction = new MeterDimerFraction(agentManager);
        meterDimerFraction.setSpeciesA(speciesA);
        meterDimerFraction.setBox(box);
        thermometer = new MeterTemperature(box, space.D());

        getController().setSleepPeriod(1);
        getController().addActivity(new ActivityIntegrate(integratorHard1));
        integratorHard1.getEventManager().addListener(new IntegratorListenerAction(new BoxImposePbc(box, space)));

        Configuration config = new ConfigurationLattice(new LatticeOrthorhombicHexagonal(space), space);
        config.initializeCoordinates(box);
    }

    public AtomLeafAgentManager<IAtom> getAgentManager() {
        return agentManager;
    }

    /**
     * Implementation of Atom.AgentSource interface. Agent is the
     * bonding partner.
     *
     * @param a ignored
     * @return Object always null
     */
    public IAtom makeAgent(IAtom a, Box agentBox) {
        return null;
    }

    public void releaseAgent(IAtom agent, IAtom atom, Box agentBox) {
    }

    public void resetBonding() {
        // nuke all bonding
        for (IAtom a : box.getLeafList()) {
            agentManager.setAgent(a, null);
        }
        IPotentialPair[][] allp = potentialMaster.getPairPotentials();
        for (int i = 0; i < box.getLeafList().size(); i++) {
            IAtom iAtom = box.getLeafList().get(i);
            int iType = iAtom.getType().getIndex();
            neighborManager.makeNeighborIterator().iterUpNeighbors(i, new NeighborIterator.NeighborConsumer() {
                @Override
                public void accept(IAtom jAtom, Vector rij) {
                    int jType = jAtom.getType().getIndex();
                    P2SquareWellBondedFasterer pij = (P2SquareWellBondedFasterer) allp[iType][jType];
                    pij.determineBonding(iAtom, jAtom, rij);

                }
            });
        }
    }
}

