/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.chainequilibrium;


import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomLeafAgentManager.AgentSource;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.config.ConfigurationLatticeRandom;
import etomica.data.meter.MeterTemperature;
import etomica.integrator.IntegratorHardFasterer;
import etomica.integrator.IntegratorMDFasterer;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.molecule.IMolecule;
import etomica.nbr.list.NeighborListManagerFastererHard;
import etomica.potential.BondingInfo;
import etomica.potential.P2HardGeneric;
import etomica.potential.compute.PotentialComputePairGeneral;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.species.SpeciesGeneral;
import etomica.units.Kelvin;

public class ChainEquilibriumFastererSim extends Simulation implements AgentSource<IAtom[]> {

    public final PotentialComputePairGeneral potentialMaster;
    public final ConfigurationLatticeRandom config;
    public IntegratorHardFasterer integratorHard;
    public java.awt.Component display;
    public Box box;
    public MeterTemperature thermometer;
    public SpeciesGeneral speciesA;
    public SpeciesGeneral speciesB;
    //    public SpeciesSpheresMono speciesC;
    public P2HardGeneric p2AA, p2BB; //, p2CC, p2BC;
    public P2SquareWellBondedFasterer ABbonded; //, ACbonded;

    public AtomLeafAgentManager<IAtom[]> agentManager;
    public int nCrossLinkersAcid;
    public int nDiol, nDiAcid;
    public int nMonoOl, nMonoAcid;

    public ChainEquilibriumFastererSim(Space space) {
        super(space);

        speciesA = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this), true);
        speciesB = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this), true);
        addSpecies(speciesA);
        addSpecies(speciesB);

        box = this.makeBox(new BoundaryRectangularPeriodic(space, space.D() == 2 ? 60 : 20));

        NeighborListManagerFastererHard neighborManager = new NeighborListManagerFastererHard(getSpeciesManager(), box, 1, 3, BondingInfo.noBonding());
        neighborManager.setDoDownNeighbors(true);
        potentialMaster = new PotentialComputePairGeneral(getSpeciesManager(), box, neighborManager);

        double diameter = 1.0;
        double lambda = 2.0;

        box.setNMolecules(speciesA, 50);
        nDiol = 50;
        box.setNMolecules(speciesB, 100);
        nDiAcid = 100;
        config = new ConfigurationLatticeRandom(space.D() == 2 ? new LatticeOrthorhombicHexagonal(space) : new LatticeCubicFcc(space), space, random);
        config.initializeCoordinates(box);

        agentManager = new AtomLeafAgentManager<>(this, box);

        //potentials
        p2AA = new P2HardGeneric(new double[]{diameter}, new double[]{Double.POSITIVE_INFINITY}, true);
        ABbonded = new P2SquareWellBondedFasterer(agentManager, diameter / lambda, lambda, 0.0);
        p2BB = new P2HardGeneric(new double[]{diameter}, new double[]{Double.POSITIVE_INFINITY}, true);

        potentialMaster.setPairPotential(speciesA.getLeafType(), speciesA.getLeafType(), p2AA);
        potentialMaster.setPairPotential(speciesA.getLeafType(), speciesB.getLeafType(), ABbonded);
        potentialMaster.setPairPotential(speciesB.getLeafType(), speciesB.getLeafType(), p2BB);

        integratorHard = new IntegratorHardFasterer(IntegratorHardFasterer.extractHardPotentials(potentialMaster), neighborManager, random, 0.002, Kelvin.UNIT.toSim(300), box);
        integratorHard.setIsothermal(true);
        integratorHard.setThermostat(IntegratorMDFasterer.ThermostatType.ANDERSEN_SINGLE);
        integratorHard.setThermostatInterval(1);

        // **** Setting Up the thermometer Meter *****

        thermometer = new MeterTemperature(box, space.D());

        getController().addActivity(new ActivityIntegrate(integratorHard, true));
    }

    public int getNMonoOl() {
        return nMonoOl;
    }

    public void setNMonoOl(int monoOl) {
        nMonoOl = monoOl;
        box.setNMolecules(speciesA, nMonoOl + nDiol);
    }

    public int getNMonoAcid() {
        return nMonoAcid;
    }

    public void setNMonoAcid(int monoAcid) {
        nMonoAcid = monoAcid;
        box.setNMolecules(speciesB, nMonoAcid + nDiAcid + nCrossLinkersAcid);
    }

    public int getNDiol() {
        return nDiol;
    }

    public void setNDiol(int diol) {
        nDiol = diol;
        box.setNMolecules(speciesA, nMonoOl + nDiol);
    }

    public int getNDiAcid() {
        return nDiAcid;
    }

    public void setNDiAcid(int diAcid) {
        nDiAcid = diAcid;
        box.setNMolecules(speciesB, nMonoAcid + nDiAcid + nCrossLinkersAcid);
    }

    public int getNCrossLinkersAcid() {
        return nCrossLinkersAcid;
    }

    public void setNCrossLinkersAcid(int crossLinkersAcid) {
        nCrossLinkersAcid = crossLinkersAcid;
        box.setNMolecules(speciesB, nMonoAcid + nDiAcid + nCrossLinkersAcid);
    }

    public void resetBonds() {
        IAtomList atoms = box.getLeafList();
        for (int i = 0; i < atoms.size(); i++) {
            IAtom a = atoms.get(i);
            agentManager.setAgent(a, makeAgent(a, box));
        }
    }

    /**
     * Implementation of AtomAgentManager.AgentSource interface. Agent
     * is used to hold bonding partners.
     */
    public IAtom[] makeAgent(IAtom a, Box agentBox) {
        IMolecule m = a.getParentGroup();
        int nBonds = 2;
        if (m.getType() == speciesA) {
            if (m.getIndex() < nMonoOl) {
                nBonds = 1;
            }
        } else {
            if (m.getIndex() < nMonoAcid) {
                nBonds = 1;
            } else if (m.getIndex() >= nMonoAcid + nDiAcid) {
                nBonds = 3;
            }
        }
        return new IAtom[nBonds];
    }

    public void releaseAgent(IAtom[] agent, IAtom atom, Box agentBox) {
    }

    public AtomLeafAgentManager<IAtom[]> getAgentManager() {
        return agentManager;
    }
}
