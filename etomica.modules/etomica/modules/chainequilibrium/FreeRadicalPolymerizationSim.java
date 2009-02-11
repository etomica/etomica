package etomica.modules.chainequilibrium;

import etomica.action.activity.ActivityIntegrate;
import etomica.api.IAtomLeaf;
import etomica.api.IAtomType;
import etomica.api.IAtomTypeSphere;
import etomica.api.IBox;
import etomica.api.IController;
import etomica.api.IMoleculeList;
import etomica.api.IPotentialMaster;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomLeafAgentManager.AgentSource;
import etomica.box.Box;
import etomica.integrator.IntegratorHard;
import etomica.integrator.IntegratorMD.ThermostatType;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.nbr.list.PotentialMasterList;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.ISpace;
import etomica.space2d.Space2D;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Kelvin;

public class FreeRadicalPolymerizationSim extends Simulation implements AgentSource {

	public IController controller1;
	public IntegratorHard integratorHard;
	public java.awt.Component display;
	public IBox box;
	public SpeciesSpheresMono speciesA; // initiator
	public SpeciesSpheresMono speciesB; // monomer
	public P2SquareWellBonded p2AA;
	public P2SquareWellRadical p2AB, p2BB;
    public ActivityIntegrate activityIntegrate;
    public AtomLeafAgentManager agentManager = null;
    public final IPotentialMaster potentialMaster;
    public final ConfigurationLatticeFreeRadical config;

    public FreeRadicalPolymerizationSim() {
        this(Space2D.getInstance());
    }
    
    public FreeRadicalPolymerizationSim(ISpace space) {
        super(space);
        potentialMaster = new PotentialMasterList(this, 3, space);
        ((PotentialMasterList)potentialMaster).setCellRange(1);

        controller1 = getController();

        double diameter = 1.0;
        double lambda = 2.0;

        integratorHard = new IntegratorHard(this, potentialMaster, space);
        integratorHard.setIsothermal(true);
        integratorHard.setTemperature(Kelvin.UNIT.toSim(300));
        integratorHard.setTimeStep(0.002);
        integratorHard.setThermostat(ThermostatType.ANDERSEN_SINGLE);
        integratorHard.setThermostatInterval(1);

        box = new Box(new BoundaryRectangularPeriodic(space, space.D() == 2 ? 60 : 20), space);
        addBox(box);
        integratorHard.setBox(box);
        integratorHard.addIntervalAction(((PotentialMasterList)potentialMaster).getNeighborManager(box));
        
        speciesA = new SpeciesSpheresMono(this, space);
        speciesB = new SpeciesSpheresMono(this, space);
        getSpeciesManager().addSpecies(speciesA);
        getSpeciesManager().addSpecies(speciesB);
        ((IAtomTypeSphere)speciesA.getLeafType()).setDiameter(diameter);
        ((IAtomTypeSphere)speciesB.getLeafType()).setDiameter(diameter);
        box.setNMolecules(speciesA, 50);
        box.setNMolecules(speciesB, 100);
        config = new ConfigurationLatticeFreeRadical(space.D() == 2 ? new LatticeOrthorhombicHexagonal(space) : new LatticeCubicFcc(space), space, random);
        config.setSpecies(speciesA, speciesB);
        config.initializeCoordinates(box);

        agentManager = new AtomLeafAgentManager(this,box);
        resetBonds();

		//potentials
        p2AA = new P2SquareWellBonded(space, agentManager, diameter / lambda, lambda, 0);
		p2AB = new P2SquareWellRadical(space, agentManager, diameter / lambda, lambda, 0.0, random);
        p2BB = new P2SquareWellRadical(space, agentManager, diameter / lambda, lambda, 0.0, random);

		potentialMaster.addPotential(p2AA,
		        new IAtomType[] { speciesA.getLeafType(), speciesA.getLeafType() });
		potentialMaster.addPotential(p2AB,
		        new IAtomType[] { speciesA.getLeafType(), speciesB.getLeafType() });
        potentialMaster.addPotential(p2BB,
                new IAtomType[] { speciesB.getLeafType(), speciesB.getLeafType() });

		// **** Setting Up the thermometer Meter *****
		
		activityIntegrate = new ActivityIntegrate(integratorHard, 1, true);
		getController().addAction(activityIntegrate);
	}
    
    public void resetBonds() {
        
        IMoleculeList initiators = box.getMoleculeList(speciesA);
        for (int i=0; i<initiators.getMoleculeCount(); i++) {
            IAtomLeaf initiator0 = initiators.getMolecule(i).getChildList().getAtom(0);
            IAtomLeaf[] bonds0 = (IAtomLeaf[])agentManager.getAgent(initiator0);
            if (i<initiators.getMoleculeCount()-1) {
                i++;
                IAtomLeaf initiator1 = initiators.getMolecule(i).getChildList().getAtom(0);
                IAtomLeaf[] bonds1 = (IAtomLeaf[])agentManager.getAgent(initiator1);
                bonds0[0] = initiator1;
                bonds1[0] = initiator0;
            }
            else {
                bonds0[0] = null;
            }
        }
        
        IMoleculeList monomers = box.getMoleculeList(speciesB);
        for (int i=0; i<monomers.getMoleculeCount(); i++) {
            IAtomLeaf[] bonds = (IAtomLeaf[])agentManager.getAgent(monomers.getMolecule(i).getChildList().getAtom(0));
            bonds[0] = null;
            bonds[1] = null;
        }
    }
    
    public Class getAgentClass() {
        return IAtomLeaf[].class;
    }
    
	/**
	 * Implementation of AtomAgentManager.AgentSource interface. Agent
     * is used to hold bonding partners.
	 */
	public Object makeAgent(IAtomLeaf a) {
	    if (a.getType() == speciesB.getLeafType()) {
	        return new IAtomLeaf[2]; // monomer
	    }
		return new IAtomLeaf[1]; // initiator
	}
    
    public void releaseAgent(Object agent, IAtomLeaf atom) {}
    
    public AtomLeafAgentManager getAgentManager() {
    	return agentManager;
    }
}
