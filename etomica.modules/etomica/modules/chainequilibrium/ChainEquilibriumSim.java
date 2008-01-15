package etomica.modules.chainequilibrium;

import etomica.action.BoxImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.action.activity.Controller;
import etomica.atom.AtomAgentManager;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomType;
import etomica.atom.AtomTypeSphere;
import etomica.atom.IAtom;
import etomica.atom.AtomAgentManager.AgentSource;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.data.meter.MeterTemperature;
import etomica.integrator.IntegratorHard;
import etomica.integrator.IntegratorMD.ThermostatType;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space2d.Space2D;
import etomica.species.Species;
import etomica.species.SpeciesSpheresMono;

public class ChainEquilibriumSim extends Simulation implements AgentSource {

	public MeterChainLength molecularCount;
	public Controller controller1;
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
    public ActivityIntegrate activityIntegrate;
    private AtomLeafAgentManager agentManager = null;
    public IAtom[] agents;

    public ChainEquilibriumSim() {
        super(Space2D.getInstance());
        PotentialMaster potentialMaster = new PotentialMaster(space);
        controller1 = getController();

        double diameter = 1.0;

        integratorHard1 = new IntegratorHard(this, potentialMaster);
        integratorHard1.setIsothermal(true);
        integratorHard1.setThermostat(ThermostatType.ANDERSEN_SINGLE);
        integratorHard1.setThermostatInterval(10);

        box = new Box(this);
        addBox(box);
        box.setBoundary(new BoundaryRectangularPeriodic(space, random, 30));
        integratorHard1.setBox(box);	
        speciesA = new SpeciesSpheresMono(this);
        speciesB = new SpeciesSpheresMono(this);
        getSpeciesManager().addSpecies(speciesA);
        getSpeciesManager().addSpecies(speciesB);
        ((AtomTypeSphere)speciesA.getLeafType()).setDiameter(diameter);
        ((AtomTypeSphere)speciesB.getLeafType()).setDiameter(diameter);
        box.setNMolecules(speciesA, 10);
        box.setNMolecules(speciesB, 40);
        new ConfigurationLattice(new LatticeOrthorhombicHexagonal()).initializeCoordinates(box);

        agentManager = new AtomLeafAgentManager(this,box);

        molecularCount = new MeterChainLength(agentManager);
        molecularCount.setBox(box);

		//potentials
        AAbonded = new P2SquareWellBonded(space, agentManager, 0.5 * diameter, 2.0, 1.0);
		ABbonded = new P2SquareWellBonded(space, agentManager, 0.5 * diameter, 2.0, 1.0);
		BBbonded = new P2SquareWellBonded(space, agentManager, 0.5 * diameter, 2.0, 1.0);

		potentialMaster.addPotential(AAbonded,
		        new AtomType[] { speciesA.getLeafType(), speciesA.getLeafType() });
		potentialMaster.addPotential(ABbonded,
		        new AtomType[] { speciesA.getLeafType(), speciesB.getLeafType() });
		
		potentialMaster.addPotential(BBbonded,
		        new AtomType[] { speciesB.getLeafType(), speciesB.getLeafType() });


		// **** Setting Up the thermometer Meter *****
		
		thermometer = new MeterTemperature();
		thermometer.setBox(box);
        
		activityIntegrate = new ActivityIntegrate(integratorHard1, 1, true);
		getController().addAction(activityIntegrate);
		integratorHard1.addIntervalAction(new BoxImposePbc(box));

	}
    
    public Class getAgentClass() {
        return IAtom[].class;
    }
    
	/**
	 * Implementation of AtomAgentManager.AgentSource interface. Agent
     * is used to hold bonding partners.
	 */
	public Object makeAgent(IAtom a) {
		
		return new IAtom[2];
	}
    
    public void releaseAgent(Object agent, IAtom atom) {}
    
    public AtomAgentManager getAgentManager() {
    	return agentManager;
    }
}
