/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.chainequilibrium;


import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomLeafAgentManager.AgentSource;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.integrator.IntegratorHardFasterer;
import etomica.integrator.IntegratorMDFasterer;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.nbr.list.NeighborListManagerFastererHard;
import etomica.potential.BondingInfo;
import etomica.potential.compute.PotentialComputePair;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.space2d.Space2D;
import etomica.species.SpeciesGeneral;
import etomica.units.Kelvin;
import etomica.util.random.RandomMersenneTwister;

public class FreeRadicalPolymerizationFastererSim extends Simulation implements AgentSource<IAtom[]> {

    public final PotentialComputePair potentialMaster;
    public final ConfigurationLatticeFreeRadical config;
    public IntegratorHardFasterer integratorHard;
    public java.awt.Component display;
    public Box box;
    public SpeciesGeneral speciesA; // initiator
    public SpeciesGeneral speciesB; // monomer
    public P2SquareWellBondedFasterer p2AA;
    public P2SquareWellRadicalFasterer p2AB, p2BB;

    public AtomLeafAgentManager<IAtom[]> agentManager = null;

    public FreeRadicalPolymerizationFastererSim() {
        this(Space2D.getInstance());
    }

    public FreeRadicalPolymerizationFastererSim(Space space) {
        super(space);
        setRandom(new RandomMersenneTwister(1));
        speciesA = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this), true);
        speciesB = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this), true);
        addSpecies(speciesA);
        addSpecies(speciesB);

        box = this.makeBox(new BoundaryRectangularPeriodic(space, space.D() == 2 ? 60 : 20));

        NeighborListManagerFastererHard neighborManager = new NeighborListManagerFastererHard(getSpeciesManager(), box, 1, 3, BondingInfo.noBonding());
        neighborManager.setDoDownNeighbors(true);
        potentialMaster = new PotentialComputePair(getSpeciesManager(), box, neighborManager);

        double diameter = 1.0;
        double lambda = 2.0;

        box.setNMolecules(speciesA, 50);
        box.setNMolecules(speciesB, 100);
        config = new ConfigurationLatticeFreeRadical(space.D() == 2 ? new LatticeOrthorhombicHexagonal(space) : new LatticeCubicFcc(space), space, random);
        config.setSpecies(speciesA, speciesB);
        config.initializeCoordinates(box);

        agentManager = new AtomLeafAgentManager<>(this, box);
        resetBonds();

        //potentials
        p2AA = new P2SquareWellBondedFasterer(agentManager, diameter / lambda, lambda, 0);
        p2AB = new P2SquareWellRadicalFasterer(agentManager, diameter / lambda, lambda, 0.0, random);
        p2BB = new P2SquareWellRadicalFasterer(agentManager, diameter / lambda, lambda, 0.0, random);

        potentialMaster.setPairPotential(speciesA.getLeafType(), speciesA.getLeafType(), p2AA);
        potentialMaster.setPairPotential(speciesA.getLeafType(), speciesB.getLeafType(), p2AB);
        potentialMaster.setPairPotential(speciesB.getLeafType(), speciesB.getLeafType(), p2BB);

        integratorHard = new IntegratorHardFasterer(IntegratorHardFasterer.extractHardPotentials(potentialMaster), neighborManager, random, 0.002, Kelvin.UNIT.toSim(300), box);
        integratorHard.setIsothermal(true);
        integratorHard.setThermostat(IntegratorMDFasterer.ThermostatType.ANDERSEN_SINGLE);
        integratorHard.setThermostatInterval(1);

        getController().addActivity(new ActivityIntegrate(integratorHard, true));
    }

    public void resetBonds() {

        IMoleculeList initiators = box.getMoleculeList(speciesA);
        for (int i = 0; i < initiators.size(); i++) {
            IAtom initiator0 = initiators.get(i).getChildList().get(0);
            IAtom[] bonds0 = agentManager.getAgent(initiator0);
            if (i < initiators.size() - 1) {
                i++;
                IAtom initiator1 = initiators.get(i).getChildList().get(0);
                IAtom[] bonds1 = agentManager.getAgent(initiator1);
                bonds0[0] = initiator1;
                bonds1[0] = initiator0;
            } else {
                bonds0[0] = null;
            }
        }

        IMoleculeList monomers = box.getMoleculeList(speciesB);
        for (IMolecule monomer : monomers) {
            IAtom[] bonds = agentManager.getAgent(monomer.getChildList().get(0));
            bonds[0] = null;
            bonds[1] = null;
        }
    }

    /**
     * Implementation of AtomAgentManager.AgentSource interface. Agent
     * is used to hold bonding partners.
     */
    public IAtom[] makeAgent(IAtom a, Box agentBox) {
        if (a.getType() == speciesB.getLeafType()) {
            return new IAtom[2]; // monomer
        }
        return new IAtom[1]; // initiator
    }

    public void releaseAgent(IAtom[] agent, IAtom atom, Box agentBox) {
    }

    public AtomLeafAgentManager<IAtom[]> getAgentManager() {
        return agentManager;
    }
}
