package etomica.modules.reactionequilibrium;

import javax.swing.JPanel;

import etomica.Atom;
import etomica.Controller;
import etomica.Default;
import etomica.Phase;
import etomica.Simulation;
import etomica.Species;
import etomica.SpeciesSpheresMono;
import etomica.action.PhaseImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.data.meter.MeterTemperature;
import etomica.graphics.DisplayPhase;
import etomica.integrator.IntegratorHard;
import etomica.space2d.Space2D;

public class ReactionEquilibrium extends Simulation implements Atom.AgentSource {

	public Controller controller1;
	public JPanel panel = new JPanel(new java.awt.BorderLayout());
	public IntegratorHard integratorHard1;
	public java.awt.Component display;
	public Phase phase1;
	public DisplayPhase displayPhase1;
	public etomica.action.SimulationRestart restartAction;
	public boolean initializing = true;
	public int idx;
	public MeterTemperature thermometer;
	public SpeciesSpheresMono speciesA;
	public SpeciesSpheresMono speciesB;
	public P2SquareWellBonded AAbonded;
	public P2SquareWellBonded ABbonded;
	public P2SquareWellBonded BBbonded;
	public MeterDimerFraction meterDimerFraction;
	public ReactionEquilibrium() {
		super(Space2D.getInstance());
        controller1 = getController();
        idx = Atom.requestAgentIndex(this);

		double diameter = 1.0;
		Default.ATOM_SIZE = diameter;

		//controller and integrator
		integratorHard1 = new IntegratorHard(potentialMaster);
		integratorHard1.setIsothermal(true);
//        integratorHard1.setThermostat(IntegratorMD.ANDERSEN_SINGLE);

		//construct phase
		phase1 = new Phase(this);
		integratorHard1.addPhase(phase1);
		speciesA = new SpeciesSpheresMono(this);
		speciesB = new SpeciesSpheresMono(this);
		speciesA.setDiameter(diameter);
        speciesA.setNMolecules(50);
        speciesB.setNMolecules(50);
        phase1.makeMolecules();

		//potentials
		AAbonded = new P2SquareWellBonded(space, idx, 0.5 * Default.ATOM_SIZE, //core
				2.0, //well multiplier
				Default.POTENTIAL_WELL);
		ABbonded = new P2SquareWellBonded(space, idx, 0.5 * Default.ATOM_SIZE, //core
				2.0, //well multiplier
				Default.POTENTIAL_WELL);
		BBbonded = new P2SquareWellBonded(space, idx, 0.5 * Default.ATOM_SIZE, //core
				2.0, //well multiplier
				Default.POTENTIAL_WELL);
/*		P2SquareWell AAbonded = new P2SquareWell(space, 0.5 * Default.ATOM_SIZE, //core
				2.0, //well multiplier
				Default.POTENTIAL_WELL);
		P2SquareWell ABbonded = new P2SquareWell(space, 0.5 * Default.ATOM_SIZE, //core
				2.0, //well multiplier
				Default.POTENTIAL_WELL);
		P2SquareWell BBbonded = new P2SquareWell(space, 0.5 * Default.ATOM_SIZE, //core
				2.0, //well multiplier
				Default.POTENTIAL_WELL);*/
		potentialMaster.setSpecies(AAbonded,
				new Species[] { speciesA, speciesA });
		potentialMaster.setSpecies(ABbonded,
				new Species[] { speciesA, speciesB });
		potentialMaster.setSpecies(BBbonded,
				new Species[] { speciesB, speciesB });

		meterDimerFraction = new MeterDimerFraction(idx);
		thermometer = new MeterTemperature();
		thermometer.setPhase(phase1);
        
		ActivityIntegrate activityIntegrate = new ActivityIntegrate(integratorHard1);
		activityIntegrate.setDoSleep(true);
		activityIntegrate.setSleepPeriod(1);
		getController().addAction(activityIntegrate);
		integratorHard1.addListener(new PhaseImposePbc(phase1));
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
		return new Atom[2];
	}

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

//=======================================================================
/**
 * A modulator specific for the square-well potential. The modulated value is
 * the ratio of the core diameter to the well diameter. This is not explicitly
 * given as a field of the potential, so the modulator must do some manipulation
 * of the accessible fields to implement changes. The value passed to setValue
 * is interpreted as 100 times this ratio, and is used to adjust the core
 * diameter and lambda value of the potential accordingly. The full well
 * diameter of the potential (lambda*sigma) is fixed at construction.
 */
