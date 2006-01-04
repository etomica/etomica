package etomica.modules.reactionequilibrium;

import javax.swing.JPanel;

import etomica.action.PhaseImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.action.activity.Controller;
import etomica.atom.Atom;
import etomica.atom.AtomAgentManager;
import etomica.atom.AtomTypeSphere;
import etomica.atom.AtomAgentManager.AgentSource;
import etomica.config.Configuration;
import etomica.config.ConfigurationSequential;
import etomica.data.meter.MeterTemperature;
import etomica.graphics.DisplayPhase;
import etomica.integrator.IntegratorHard;
import etomica.phase.Phase;
import etomica.simulation.Simulation;
import etomica.space2d.Space2D;
import etomica.space3d.Space3D;
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
    public AtomAgentManager agentManager;
    public Atom[] agents;
    
    public ReactionEquilibrium() {
        super(Space2D.getInstance());
        defaults.ignoreOverlap = true;
        controller1 = getController();

        double diameter = 1.0;
        defaults.atomSize = diameter;

        //controller and integrator
        integratorHard1 = new IntegratorHard(this);
        integratorHard1.setIsothermal(true);
//        integratorHard1.setThermostat(IntegratorMD.ANDERSEN_SINGLE);

        //construct phase
        phase1 = new Phase(this);
        integratorHard1.setPhase(phase1);
        speciesA = new SpeciesSpheresMono(this);
        speciesB = new SpeciesSpheresMono(this);
        ((AtomTypeSphere)speciesA.getMoleculeType()).setDiameter(diameter);
        speciesA.setNMolecules(30);
        speciesB.setNMolecules(30);
        phase1.makeMolecules();
        Configuration config = new ConfigurationSequential(space);
        config.initializeCoordinates(phase1);

        //potentials
        AAbonded = new P2SquareWellBonded(space, this, 0.5 * defaults.atomSize, //core
                2.0, //well multiplier
                defaults.potentialWell, defaults.ignoreOverlap);
        ABbonded = new P2SquareWellBonded(space, this, 0.5 * defaults.atomSize, //core
                2.0, //well multiplier
                defaults.potentialWell, defaults.ignoreOverlap);
        BBbonded = new P2SquareWellBonded(space, this, 0.5 * defaults.atomSize, //core
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

        meterDimerFraction = new MeterDimerFraction(this);
        meterDimerFraction.setPhase(phase1);
        thermometer = new MeterTemperature();
        thermometer.setPhase(phase1);
        
        ActivityIntegrate activityIntegrate = new ActivityIntegrate(this,integratorHard1);
        activityIntegrate.setDoSleep(true);
        activityIntegrate.setSleepPeriod(1);
        getController().addAction(activityIntegrate);
        integratorHard1.addListener(new PhaseImposePbc(phase1));
	}
    
    public Atom[] getAgents(Phase phase) {
        // the other classes don't know it, but there's only one phase.  :)
        if (agentManager == null) {
          agentManager = new AtomAgentManager(this,phase);
        }
        return (Atom[])agentManager.getAgents();
    }

    /**
     * Implementation of Atom.AgentSource interface, returning null. Agent in
     * atom is used to hold bonding partner.
     * 
     * @param a
     *            ignored
     * @return Object always null
     */
    public Object makeAgent(Atom a) {
        if (a == null) {
            // return bogus atom so AtomAgentManager knows the class
            return new Atom(Space3D.getInstance());
        }
        return null;
    }
    
    public void releaseAgent(Object agent) {}

	public static void main(String[] args) {
		javax.swing.JFrame f = new javax.swing.JFrame(); //create a window
		f.setSize(800, 550);

		ReactionEquilibrium sim = new ReactionEquilibrium();
		ReactionEquilibriumGraphic graphic = new ReactionEquilibriumGraphic(sim);
		f.getContentPane().add(graphic.panel);
		f.pack();
		f.show();
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

}//end of KineticsModule class

