package etomica.modules.chainequilibrium;

import javax.swing.JPanel;

import etomica.action.PhaseImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.action.activity.Controller;
import etomica.atom.Atom;
import etomica.config.ConfigurationSequential;
import etomica.data.meter.MeterTemperature;
import etomica.integrator.IntegratorHard;
import etomica.phase.Phase;
import etomica.simulation.Simulation;
import etomica.space2d.Space2D;
import etomica.species.Species;
import etomica.species.SpeciesSpheresMono;

public class ReactionEquilibrium extends Simulation implements Atom.AgentSource {

	public MeterChainLength molecularCount;
	public Controller controller1;
	public JPanel panel = new JPanel(new java.awt.BorderLayout());
	public IntegratorHard integratorHard1;
	public java.awt.Component display;
	public Phase phase1;
	public etomica.action.SimulationRestart restartAction;
	public boolean initializing = true;
	public int idx;
	public MeterTemperature thermometer;
	public SpeciesSpheresMono speciesA;
	public SpeciesSpheresMono speciesB;
	public P2SquareWellBonded AAbonded;
	public P2SquareWellBonded ABbonded;
	public P2SquareWellBonded BBbonded;
	
    public ReactionEquilibrium() {
        super(Space2D.getInstance());
        controller1 = getController();

        double diameter = 1.0;
        idx = Atom.requestAgentIndex(this);
        
        molecularCount = new MeterChainLength(idx);
		

        getDefaults().atomSize = diameter;

        integratorHard1 = new IntegratorHard(this);
        integratorHard1.setIsothermal(true);

        phase1 = new Phase(this);
        integratorHard1.addPhase(phase1);	
        speciesA = new SpeciesSpheresMono(this);
        speciesB = new SpeciesSpheresMono(this);
        speciesA.setDiameter(diameter);
        speciesA.setNMolecules(10);      
        speciesB.setNMolecules(40);      
        phase1.makeMolecules();
        new ConfigurationSequential(space).initializeCoordinates(phase1);

        molecularCount.setPhase(phase1);
    	
		//potentials
        AAbonded = new P2SquareWellBonded(space, idx, 0.5 * getDefaults().atomSize, 
                2.0, getDefaults().potentialWell);
		ABbonded = new P2SquareWellBonded(space, idx, 0.5 * getDefaults().atomSize,
		        2.0, getDefaults().potentialWell);
		BBbonded = new P2SquareWellBonded(space, idx, 0.5 * getDefaults().atomSize,
		        2.0, getDefaults().potentialWell);

		potentialMaster.setSpecies(AAbonded,
		        new Species[] { speciesA, speciesA });
		potentialMaster.setSpecies(ABbonded,
		        new Species[] { speciesA, speciesB });
		
		potentialMaster.setSpecies(BBbonded,new Species[] { speciesB, speciesB });


		// **** Setting Up the thermometer Meter *****
		
		thermometer = new MeterTemperature();
		thermometer.setPhase(phase1);
        
		ActivityIntegrate activityIntegrate = new ActivityIntegrate(integratorHard1,true,true);
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
		
		// ******* MARKER ******
//		System.out.println("makeAgent: 4 atom array in reactionequlibrum");
		
		return new Atom[4];
	}
	public static void main(String[] args) {
		javax.swing.JFrame f = new javax.swing.JFrame(); //create a window
		f.setSize(800, 550);
		
		// ******* MARKER ******
		System.out.println("New reaction Equilibrium");
		ReactionEquilibrium sim = new ReactionEquilibrium();
		
		// ******* MARKER ******
		System.out.println("New Reaction equilbrium graphic");
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
