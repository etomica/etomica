/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.chainequilibrium;

import etomica.action.activity.ActivityIntegrate;
import etomica.action.activity.Controller;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomLeafAgentManager.AgentSource;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.meter.MeterTemperature;
import etomica.integrator.IntegratorHard;
import etomica.integrator.IntegratorMD.ThermostatType;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.molecule.IMolecule;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.P2HardSphere;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Kelvin;

public class ChainEquilibriumSim extends Simulation implements AgentSource<IAtom[]> {

    public final PotentialMaster potentialMaster;
    public final ConfigurationLatticeRandom config;
    public Controller controller1;
	public IntegratorHard integratorHard;
	public java.awt.Component display;
	public Box box;
	public MeterTemperature thermometer;
	public SpeciesSpheresMono speciesA;
	public SpeciesSpheresMono speciesB;
//    public SpeciesSpheresMono speciesC;
	public P2HardSphere p2AA, p2BB; //, p2CC, p2BC;
	public P2SquareWellBonded ABbonded; //, ACbonded;
    public ActivityIntegrate activityIntegrate;
    public AtomLeafAgentManager<IAtom[]> agentManager = null;
    public int nCrossLinkersAcid;
    public int nDiol, nDiAcid;
    public int nMonoOl, nMonoAcid;

    public ChainEquilibriumSim(Space space) {
        super(space);

        speciesA = new SpeciesSpheresMono(this, space);
        speciesA.setIsDynamic(true);
        speciesB = new SpeciesSpheresMono(this, space);
        speciesB.setIsDynamic(true);
        addSpecies(speciesA);
        addSpecies(speciesB);

        potentialMaster = new PotentialMasterList(this, 3, space);
        ((PotentialMasterList) potentialMaster).setCellRange(1);

        controller1 = getController();

        double diameter = 1.0;
        double lambda = 2.0;


        box = this.makeBox(new BoundaryRectangularPeriodic(space, space.D() == 2 ? 60 : 20));
        box.setNMolecules(speciesA, 50);
        nDiol = 50;
        box.setNMolecules(speciesB, 100);
        nDiAcid = 100;
        config = new ConfigurationLatticeRandom(space.D() == 2 ? new LatticeOrthorhombicHexagonal(space) : new LatticeCubicFcc(space), space, random);
        config.initializeCoordinates(box);

        integratorHard = new IntegratorHard(this, potentialMaster, box);
        integratorHard.setIsothermal(true);
        integratorHard.setTemperature(Kelvin.UNIT.toSim(300));
        integratorHard.setTimeStep(0.002);
        integratorHard.setThermostat(ThermostatType.ANDERSEN_SINGLE);
        integratorHard.setThermostatInterval(1);
        integratorHard.getEventManager().addListener(((PotentialMasterList) potentialMaster).getNeighborManager(box));

        agentManager = new AtomLeafAgentManager<IAtom[]>(this, box);

        //potentials
        p2AA = new P2HardSphere(space, diameter, true);
        ABbonded = new P2SquareWellBonded(space, agentManager, diameter / lambda, lambda, 0.0);
        p2BB = new P2HardSphere(space, diameter, true);

        potentialMaster.addPotential(p2AA,
                new AtomType[]{speciesA.getLeafType(), speciesA.getLeafType()});
        potentialMaster.addPotential(ABbonded,
                new AtomType[]{speciesA.getLeafType(), speciesB.getLeafType()});

        potentialMaster.addPotential(p2BB,
                new AtomType[]{speciesB.getLeafType(), speciesB.getLeafType()});

        // **** Setting Up the thermometer Meter *****

        thermometer = new MeterTemperature(box, space.D());

        activityIntegrate = new ActivityIntegrate(integratorHard, 1, true);
        getController().addAction(activityIntegrate);
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
        for (int i = 0; i<atoms.size(); i++) {
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
	    }
	    else {
	        if (m.getIndex() < nMonoAcid) {
	            nBonds = 1;
	        }
	        else if (m.getIndex() >= nMonoAcid+nDiAcid) {
	            nBonds = 3;
	        }
	    }
		return new IAtom[nBonds];
	}
    
    public void releaseAgent(IAtom[] agent, IAtom atom, Box agentBox) {}
    
    public AtomLeafAgentManager<IAtom[]> getAgentManager() {
    	return agentManager;
    }
}
