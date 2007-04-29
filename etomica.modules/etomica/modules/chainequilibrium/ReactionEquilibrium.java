package etomica.modules.chainequilibrium;

import javax.swing.JPanel;

import etomica.action.PhaseImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.action.activity.Controller;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomTypeSphere;
import etomica.atom.IAtom;
import etomica.atom.AtomAgentManager.AgentSource;
import etomica.config.ConfigurationLattice;
import etomica.data.meter.MeterTemperature;
import etomica.integrator.IntegratorHard;
import etomica.integrator.IntervalActionAdapter;
import etomica.integrator.IntegratorMD.ThermostatType;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.phase.Phase;
import etomica.simulation.Simulation;
import etomica.space2d.Space2D;
import etomica.species.Species;
import etomica.species.SpeciesSpheresMono;

public class ReactionEquilibrium extends Simulation implements AgentSource {

	public MeterChainLength molecularCount;
	public Controller controller1;
	public JPanel panel = new JPanel(new java.awt.BorderLayout());
	public IntegratorHard integratorHard1;
	public java.awt.Component display;
	public Phase phase1;
	public etomica.action.SimulationRestart restartAction;
	public boolean initializing = true;
	public MeterTemperature thermometer;
	public SpeciesSpheresMono speciesA;
	public SpeciesSpheresMono speciesB;
	public P2SquareWellBonded AAbonded;
	public P2SquareWellBonded ABbonded;
	public P2SquareWellBonded BBbonded;
    public AtomLeafAgentManager agentManager;
    public IAtom[] agents;
	
    public ReactionEquilibrium() {
        super(Space2D.getInstance());
        controller1 = getController();

        double diameter = 1.0;
        
        molecularCount = new MeterChainLength(this);
		

        getDefaults().atomSize = diameter;

        integratorHard1 = new IntegratorHard(this);
        integratorHard1.setIsothermal(true);
        integratorHard1.setThermostat(ThermostatType.ANDERSEN_SINGLE);
        integratorHard1.setThermostatInterval(10);

        phase1 = new Phase(this);
        integratorHard1.setPhase(phase1);	
        speciesA = new SpeciesSpheresMono(this);
        speciesB = new SpeciesSpheresMono(this);
        getSpeciesManager().addSpecies(speciesA);
        getSpeciesManager().addSpecies(speciesB);
        ((AtomTypeSphere)speciesA.getMoleculeType()).setDiameter(diameter);
        phase1.getAgent(speciesA).setNMolecules(10);
        phase1.getAgent(speciesB).setNMolecules(40);
        new ConfigurationLattice(new LatticeOrthorhombicHexagonal()).initializeCoordinates(phase1);

        molecularCount.setPhase(phase1);
    	
		//potentials
        AAbonded = new P2SquareWellBonded(space, this, 0.5 * getDefaults().atomSize, 
                2.0, getDefaults().potentialWell);
		ABbonded = new P2SquareWellBonded(space, this, 0.5 * getDefaults().atomSize,
		        2.0, getDefaults().potentialWell);
		BBbonded = new P2SquareWellBonded(space, this, 0.5 * getDefaults().atomSize,
		        2.0, getDefaults().potentialWell);

		potentialMaster.addPotential(AAbonded,
		        new Species[] { speciesA, speciesA });
		potentialMaster.addPotential(ABbonded,
		        new Species[] { speciesA, speciesB });
		
		potentialMaster.addPotential(BBbonded,new Species[] { speciesB, speciesB });


		// **** Setting Up the thermometer Meter *****
		
		thermometer = new MeterTemperature();
		thermometer.setPhase(phase1);
        
		ActivityIntegrate activityIntegrate = new ActivityIntegrate(integratorHard1,true,true);
		activityIntegrate.setDoSleep(true);
		activityIntegrate.setSleepPeriod(1);
		getController().addAction(activityIntegrate);
		integratorHard1.addListener(new IntervalActionAdapter(new PhaseImposePbc(phase1)));
        agentManager = new AtomLeafAgentManager(this,phase1);
	}
    
    public AtomLeafAgentManager getAgentManager() {
        return agentManager;
    }

    public Class getAgentClass() {
        return IAtom[].class;
    }
    
	/**
	 * Implementation of AtomAgentManager.AgentSource interface. Agent
     * is used to hold bonding partners.
	 */
	public Object makeAgent(IAtom a) {
		
		return new IAtom[4];
	}
    
    public void releaseAgent(Object agent, IAtom atom) {}
}
