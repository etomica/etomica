package etomica.modules.reactionequilibrium;

import javax.swing.JFrame;
import javax.swing.JPanel;

import etomica.action.PhaseImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.action.activity.Controller;
import etomica.atom.AtomAgentManager;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomTypeSphere;
import etomica.atom.IAtom;
import etomica.atom.AtomAgentManager.AgentSource;
import etomica.config.Configuration;
import etomica.config.ConfigurationLattice;
import etomica.data.meter.MeterTemperature;
import etomica.graphics.DisplayPhase;
import etomica.integrator.IntegratorHard;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.phase.Phase;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space2d.Space2D;
import etomica.species.Species;
import etomica.species.SpeciesSpheresMono;

public class ReactionEquilibrium extends Simulation implements AgentSource {

    public Controller controller1;
    public JPanel panel = new JPanel(new java.awt.BorderLayout());
    public IntegratorHard integratorHard1;
    public java.awt.Component display;
    public Phase phase1;
    public DisplayPhase displayPhase1;
    public etomica.action.SimulationRestart restartAction;
    public boolean initializing = true;
    public MeterTemperature thermometer;
    public SpeciesSpheresMono speciesA;
    public SpeciesSpheresMono speciesB;
    public P2SquareWellBonded AAbonded;
    public P2SquareWellBonded ABbonded;
    public P2SquareWellBonded BBbonded;
    public MeterDimerFraction meterDimerFraction;
    private AtomLeafAgentManager agentManager = null;
    public IAtom[] agents;
    
    public ReactionEquilibrium() {
        super(Space2D.getInstance());
        PotentialMaster potentialMaster = new PotentialMaster(this);
        defaults.ignoreOverlap = true;
        controller1 = getController();

        double diameter = 1.0;
        defaults.atomSize = diameter;

        //controller and integrator
        integratorHard1 = new IntegratorHard(this, potentialMaster);
        integratorHard1.setIsothermal(true);
//        integratorHard1.setThermostat(IntegratorMD.ANDERSEN_SINGLE);

        //construct phase
        phase1 = new Phase(this);
        addPhase(phase1);
        integratorHard1.setPhase(phase1);
        speciesA = new SpeciesSpheresMono(this);
        speciesB = new SpeciesSpheresMono(this);
        getSpeciesManager().addSpecies(speciesA);
        getSpeciesManager().addSpecies(speciesB);
        ((AtomTypeSphere)speciesA.getMoleculeType()).setDiameter(diameter);
        phase1.getAgent(speciesA).setNMolecules(30);
        phase1.getAgent(speciesB).setNMolecules(30);
        Configuration config = new ConfigurationLattice(new LatticeOrthorhombicHexagonal());
        config.initializeCoordinates(phase1);

        agentManager = new AtomLeafAgentManager(this,phase1);

        //potentials
        AAbonded = new P2SquareWellBonded(space, agentManager, 0.5 * defaults.atomSize, //core
                2.0, //well multiplier
                defaults.potentialWell, defaults.ignoreOverlap);
        ABbonded = new P2SquareWellBonded(space, agentManager, 0.5 * defaults.atomSize, //core
                2.0, //well multiplier
                defaults.potentialWell, defaults.ignoreOverlap);
        BBbonded = new P2SquareWellBonded(space, agentManager, 0.5 * defaults.atomSize, //core
                2.0, //well multiplier
                defaults.potentialWell, defaults.ignoreOverlap);
/*      P2SquareWell AAbonded = new P2SquareWell(space, 0.5 * Default.atomSize, //core
                2.0, //well multiplier
                Default.POTENTIAL_WELL);
        P2SquareWell ABbonded = new P2SquareWell(space, 0.5 * Default.atomSize, //core
                2.0, //well multiplier
                Default.POTENTIAL_WELL);
        P2SquareWell BBbonded = new P2SquareWell(space, 0.5 * Default.atomSize, //core
                2.0, //well multiplier
                Default.POTENTIAL_WELL);*/
        potentialMaster.addPotential(AAbonded,
                new Species[] { speciesA, speciesA });
        potentialMaster.addPotential(ABbonded,
                new Species[] { speciesA, speciesB });
        potentialMaster.addPotential(BBbonded,
                new Species[] { speciesB, speciesB });

        meterDimerFraction = new MeterDimerFraction(agentManager);
        meterDimerFraction.setSpeciesA(speciesA);
        meterDimerFraction.setPhase(phase1);
        thermometer = new MeterTemperature();
        thermometer.setPhase(phase1);
        
        ActivityIntegrate activityIntegrate = new ActivityIntegrate(this,integratorHard1);
        activityIntegrate.setDoSleep(true);
        activityIntegrate.setSleepPeriod(1);
        getController().addAction(activityIntegrate);
        integratorHard1.addIntervalAction(new PhaseImposePbc(phase1));

	}
    
    public Class getAgentClass() {
        return IAtom.class;
    }

    public AtomAgentManager getAgentManager() {
    	return agentManager;
    }

    /**
     * Implementation of Atom.AgentSource interface. Agent is the 
     * bonding partner.
     * 
     * @param a  ignored
     * @return Object always null
     */
    public Object makeAgent(IAtom a) {
        return null;
    }
    
    public void releaseAgent(Object agent, IAtom atom) {}

	public static void main(String[] args) {
		JFrame f = new JFrame("Reaction Equilibrium"); //create a window
		f.setSize(800, 550);

		ReactionEquilibrium sim = new ReactionEquilibrium();
		ReactionEquilibriumGraphic graphic = new ReactionEquilibriumGraphic(sim);
		f.getContentPane().add(graphic.graphic());
		f.pack();
		f.setVisible(true);
		f.addWindowListener(new java.awt.event.WindowAdapter() { //anonymous
					// class to
					// handle
					// window
					// closing
					public void windowClosing(java.awt.event.WindowEvent e) {
						System.exit(0);
					}
				});
		//     sim.controller1.start();
	}//end of main

}

