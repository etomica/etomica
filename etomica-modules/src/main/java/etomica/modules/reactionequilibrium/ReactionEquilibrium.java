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
import etomica.integrator.IntegratorHard;
import etomica.integrator.IntegratorListenerAction;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.potential.P1HardPeriodic;
import etomica.potential.PotentialMaster;
import etomica.potential.PotentialMasterMonatomic;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space2d.Space2D;
import etomica.species.SpeciesSpheresMono;

import javax.swing.*;

public class ReactionEquilibrium extends Simulation implements AgentSource<IAtom> {

    public JPanel panel = new JPanel(new java.awt.BorderLayout());
    public IntegratorHard integratorHard1;
    public java.awt.Component display;
    public Box box;
    public etomica.action.SimulationRestart restartAction;
    public boolean initializing = true;
    public MeterTemperature thermometer;
    public SpeciesSpheresMono speciesA;
    public SpeciesSpheresMono speciesB;
    public P2SquareWellBonded AAbonded;
    public P2SquareWellBonded ABbonded;
    public P2SquareWellBonded BBbonded;
    public MeterDimerFraction meterDimerFraction;

    public IAtom[] agents;
    private AtomLeafAgentManager<IAtom> agentManager = null;
    
    public ReactionEquilibrium() {
        super(Space2D.getInstance());

        speciesA = new SpeciesSpheresMono(this, space);
        speciesA.setIsDynamic(true);
        speciesB = new SpeciesSpheresMono(this, space);
        speciesB.setIsDynamic(true);
        addSpecies(speciesA);
        addSpecies(speciesB);

        PotentialMaster potentialMaster = new PotentialMasterMonatomic(this);

        double diameter = 1.0;

        //controller and integrator
        box = this.makeBox(new BoundaryRectangularPeriodic(space, 30.0));
        integratorHard1 = new IntegratorHard(this, potentialMaster, box);
        integratorHard1.setIsothermal(true);

        box.setNMolecules(speciesA, 30);
        box.setNMolecules(speciesB, 30);
        P1HardPeriodic nullPotential = new P1HardPeriodic(space, diameter);
        integratorHard1.setNullPotential(nullPotential, speciesA.getLeafType());
        integratorHard1.setNullPotential(nullPotential, speciesB.getLeafType());

        agentManager = new AtomLeafAgentManager<IAtom>(this, box);

        //potentials
        AAbonded = new P2SquareWellBonded(space, agentManager, 0.5 * diameter, //core
                2.0, //well multiplier
                1.0, true);
        ABbonded = new P2SquareWellBonded(space, agentManager, 0.5 * diameter, //core
                2.0, //well multiplier
                1.0, true);
        BBbonded = new P2SquareWellBonded(space, agentManager, 0.5 * diameter, //core
                2.0, //well multiplier
                1.0, true);
        potentialMaster.addPotential(AAbonded,
                new AtomType[]{speciesA.getLeafType(), speciesA.getLeafType()});
        potentialMaster.addPotential(ABbonded,
                new AtomType[]{speciesA.getLeafType(), speciesB.getLeafType()});
        potentialMaster.addPotential(BBbonded,
                new AtomType[]{speciesB.getLeafType(), speciesB.getLeafType()});

        meterDimerFraction = new MeterDimerFraction(agentManager);
        meterDimerFraction.setSpeciesA(speciesA);
        meterDimerFraction.setBox(box);
        thermometer = new MeterTemperature(box, space.D());

        getController().addActivity(new ActivityIntegrate(integratorHard1)).setSleepPeriod(1);
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
     * @param a  ignored
     * @return Object always null
     */
    public IAtom makeAgent(IAtom a, Box agentBox) {
        return null;
    }
    
    public void releaseAgent(IAtom agent, IAtom atom, Box agentBox) {}


}

