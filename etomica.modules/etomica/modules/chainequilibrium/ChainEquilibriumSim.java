package etomica.modules.chainequilibrium;

import etomica.action.BoxImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.api.IAtomLeaf;
import etomica.api.IAtomTypeLeaf;
import etomica.api.IAtomTypeSphere;
import etomica.api.IBox;
import etomica.api.IController;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomLeafAgentManager.AgentSource;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.data.meter.MeterTemperature;
import etomica.integrator.IntegratorHard;
import etomica.integrator.IntegratorMD.ThermostatType;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.potential.P1HardPeriodic;
import etomica.potential.P2HardSphere;
import etomica.potential.PotentialMasterMonatomic;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space2d.Space2D;
import etomica.species.SpeciesSpheresMono;

public class ChainEquilibriumSim extends Simulation implements AgentSource {

	public IController controller1;
	public IntegratorHard integratorHard1;
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

    public ChainEquilibriumSim() {
        super(Space2D.getInstance());
        PotentialMasterMonatomic potentialMaster = new PotentialMasterMonatomic(this, space);
        controller1 = getController();

        double diameter = 1.0;

        integratorHard1 = new IntegratorHard(this, potentialMaster, space);
        integratorHard1.setIsothermal(true);
        integratorHard1.setThermostat(ThermostatType.ANDERSEN);
        integratorHard1.setThermostatInterval(100);

        box = new Box(this, space);
        addBox(box);
        box.setBoundary(new BoundaryRectangularPeriodic(space, random, 30));
        integratorHard1.setBox(box);	
        speciesA = new SpeciesSpheresMono(this, space);
        speciesB = new SpeciesSpheresMono(this, space);
        speciesC = new SpeciesSpheresMono(this, space);
        getSpeciesManager().addSpecies(speciesA);
        getSpeciesManager().addSpecies(speciesB);
        getSpeciesManager().addSpecies(speciesC);
        ((IAtomTypeSphere)speciesA.getLeafType()).setDiameter(diameter);
        ((IAtomTypeSphere)speciesB.getLeafType()).setDiameter(diameter);
        ((IAtomTypeSphere)speciesC.getLeafType()).setDiameter(diameter);
        box.setNMolecules(speciesA, 10);
        box.setNMolecules(speciesB, 40);
        box.setNMolecules(speciesC, 0);
        new ConfigurationLattice(new LatticeOrthorhombicHexagonal(), space).initializeCoordinates(box);

        agentManager = new AtomLeafAgentManager(this,box);

		//potentials
        p2AA = new P2HardSphere(space, diameter, true);
		ABbonded = new P2SquareWellBonded(space, agentManager, 0.5 * diameter, 2.0, 1.0);
        p2BB = new P2HardSphere(space, diameter, true);
        ACbonded = new P2SquareWellBonded(space, agentManager, 0.5 * diameter, 2.0, 1.0);
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
		
		integratorHard1.setNullPotential(new P1HardPeriodic(space, 3), speciesA.getLeafType());
        integratorHard1.setNullPotential(new P1HardPeriodic(space, 3), speciesB.getLeafType());
        
		activityIntegrate = new ActivityIntegrate(integratorHard1, 1, true);
		getController().addAction(activityIntegrate);
		integratorHard1.addIntervalAction(new BoxImposePbc(box, space));

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
