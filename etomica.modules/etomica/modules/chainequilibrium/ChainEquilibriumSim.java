package etomica.modules.chainequilibrium;

import etomica.action.activity.ActivityIntegrate;
import etomica.api.IAtomLeaf;
import etomica.api.IAtomTypeLeaf;
import etomica.api.IAtomTypeSphere;
import etomica.api.IBox;
import etomica.api.IController;
import etomica.api.IPotentialMaster;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomLeafAgentManager.AgentSource;
import etomica.box.Box;
import etomica.data.meter.MeterTemperature;
import etomica.integrator.IntegratorHard;
import etomica.integrator.IntegratorMD.ThermostatType;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.P2HardSphere;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space2d.Space2D;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Kelvin;

public class ChainEquilibriumSim extends Simulation implements AgentSource {

	public IController controller1;
	public IntegratorHard integratorHard;
	public java.awt.Component display;
	public IBox box;
	public MeterTemperature thermometer;
	public SpeciesSpheresMono speciesA;
	public SpeciesSpheresMono speciesB;
    public SpeciesSpheresMono speciesC;
	public P2HardSphere p2AA, p2BB, p2CC, p2BC;
	public P2SquareWellBonded ABbonded, ACbonded;
    public ActivityIntegrate activityIntegrate;
    public AtomLeafAgentManager agentManager = null;
    public final IPotentialMaster potentialMaster;
    public final ConfigurationLatticeRandom config;

    public ChainEquilibriumSim() {
        super(Space2D.getInstance());
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

        box = new Box(new BoundaryRectangularPeriodic(space, random, 60), space);
        addBox(box);
        integratorHard.setBox(box);
        integratorHard.addNonintervalListener(((PotentialMasterList)potentialMaster).getNeighborManager(box));
        integratorHard.addIntervalAction(((PotentialMasterList)potentialMaster).getNeighborManager(box));
        
        speciesA = new SpeciesSpheresMono(this, space);
        speciesB = new SpeciesSpheresMono(this, space);
        speciesC = new SpeciesSpheresMono(this, space);
        getSpeciesManager().addSpecies(speciesA);
        getSpeciesManager().addSpecies(speciesB);
        getSpeciesManager().addSpecies(speciesC);
        ((IAtomTypeSphere)speciesA.getLeafType()).setDiameter(diameter);
        ((IAtomTypeSphere)speciesB.getLeafType()).setDiameter(diameter);
        ((IAtomTypeSphere)speciesC.getLeafType()).setDiameter(diameter);
        box.setNMolecules(speciesA, 50);
        box.setNMolecules(speciesB, 100);
        box.setNMolecules(speciesC, 0);
        config = new ConfigurationLatticeRandom(new LatticeOrthorhombicHexagonal(), space, random);
        config.initializeCoordinates(box);

        agentManager = new AtomLeafAgentManager(this,box);

		//potentials
        p2AA = new P2HardSphere(space, diameter, true);
		ABbonded = new P2SquareWellBonded(space, agentManager, diameter / lambda, lambda, 0.0);
        p2BB = new P2HardSphere(space, diameter, true);
        ACbonded = new P2SquareWellBonded(space, agentManager, diameter / lambda, lambda, 0.0);
        p2BC = new P2HardSphere(space, diameter, true);
        p2CC = new P2HardSphere(space, diameter, true);

		potentialMaster.addPotential(p2AA,
		        new IAtomTypeLeaf[] { speciesA.getLeafType(), speciesA.getLeafType() });
		potentialMaster.addPotential(ABbonded,
		        new IAtomTypeLeaf[] { speciesA.getLeafType(), speciesB.getLeafType() });
        potentialMaster.addPotential(ACbonded,
                new IAtomTypeLeaf[] { speciesA.getLeafType(), speciesC.getLeafType() });
		
		potentialMaster.addPotential(p2BB,
		        new IAtomTypeLeaf[] { speciesB.getLeafType(), speciesB.getLeafType() });
        potentialMaster.addPotential(p2BC,
                new IAtomTypeLeaf[] { speciesB.getLeafType(), speciesC.getLeafType() });

        potentialMaster.addPotential(p2CC,
                new IAtomTypeLeaf[] { speciesC.getLeafType(), speciesC.getLeafType() });

		// **** Setting Up the thermometer Meter *****
		
		thermometer = new MeterTemperature(box, space.D());

		activityIntegrate = new ActivityIntegrate(integratorHard, 1, true);
		getController().addAction(activityIntegrate);
	}
    
    public Class getAgentClass() {
        return IAtomLeaf[].class;
    }
    
	/**
	 * Implementation of AtomAgentManager.AgentSource interface. Agent
     * is used to hold bonding partners.
	 */
	public Object makeAgent(IAtomLeaf a) {
	    if (a.getType() == speciesC.getLeafType()) {
	        return new IAtomLeaf[3];
	    }
		return new IAtomLeaf[2];
	}
    
    public void releaseAgent(Object agent, IAtomLeaf atom) {}
    
    public AtomLeafAgentManager getAgentManager() {
    	return agentManager;
    }
}
