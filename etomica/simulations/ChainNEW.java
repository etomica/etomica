package etomica.simulations;
import java.awt.Component;

import etomica.*;
import etomica.graphics.*;
import java.awt.*;
import javax.swing.*;

/**
 * Simple hard-sphere molecular dynamics simulation in 2D.
 *
 * @author David Kofke
 */
 
public class ChainNEW extends SimulationGraphic {
    
	public IntegratorHard integrator;
	public SpeciesSpheresMono species;
	public Phase phase;
	public Potential2 potential;
	public Controller controller;
	public DisplayPhase display;

	public ChainNEW() {
  //      super(new etomica.space.continuum.Space(2));
		super(new Space2D());
		Default.ATOM_SIZE = 1.0;
 //can't use cell list until integrator is updated for it      setIteratorFactory(new IteratorFactoryCell(this));
		Simulation.instance = this;
		integrator = new IntegratorHard(this);
//		integrator.setIsothermal(true);
		species = new SpeciesSpheresMono(this);
		SpeciesSpheres speciesChain = new SpeciesSpheres(1,10);
		species.setNMolecules(100);
		phase = new Phase(this);
		potential = new P2HardSphere();
		potential.setSpecies(species,species);
		P1TetheredHardSpheres potentialChainIntra = new P1TetheredHardSpheres();
		PotentialGroup p2Inter = new PotentialGroup(2);
		P2HardSphere chainSphere = new P2HardSphere(p2Inter);
		potentialChainIntra.p2Tether.setTetherLength(Default.ATOM_SIZE);
		potentialChainIntra.setSpecies(new Species[] {speciesChain});
		p2Inter.setSpecies(new Species[] {species, speciesChain});
		controller = new Controller(this);
		new DeviceTrioControllerButton(this, controller);
		display = new DisplayPhase(this);
		DisplayTimer timer = new DisplayTimer(integrator);
		timer.setUpdateInterval(10);
		ColorSchemeByType.setColor(species, java.awt.Color.red);
		ColorSchemeByType.setColor(((AtomFactoryMono)((AtomFactoryHomo)speciesChain.moleculeFactory()).childFactory()).type(),java.awt.Color.blue);
//		panel().setBackground(java.awt.Color.yellow);
		DeviceKicker dk = new DeviceKicker(display);
		elementCoordinator.go();
		
		Atom first = phase.getAgent(speciesChain).node.firstLeafAtom();
		Atom last = phase.getAgent(speciesChain).node.lastLeafAtom();
		first.coord.momentum().E(0.0);
		first.coord.setMass(Double.MAX_VALUE);
		Atom chain = ((AtomTreeNodeGroup)phase.getAgent(speciesChain).node).firstChildAtom();
		chain.coord.translateBy(new Space2D.Vector(4,15));
		dk.setAtom(last);
		
	}
    
	/**
	 * Demonstrates how this class is implemented.
	 */
	public static void main(String[] args) {
		ChainNEW sim = new ChainNEW();
		SimulationGraphic.makeAndDisplayFrame(sim);
	//	sim.controller.start();
	}//end of main
	
	

	public static class DeviceKicker extends Device implements Drawable{
    
		private DisplayPhase display;
		private boolean visible=false;
		private Atom target;
		private Space.Vector impulse;
		private Space2D.Vector targetpos;
		private Space2D.Vector mousepos;
		private double impulsefactor;
		private double impulsefactor1 = .0051;

		public DeviceKicker(DisplayPhase display) {
			this(Simulation.instance, display);
		}
		public DeviceKicker(Simulation sim, DisplayPhase display) {
			super(sim);
			this.display = display;
			display.addDisplayPhaseListener(new Kicker());
			display.addDrawable(this);
			impulse = sim.space.makeVector();
			mousepos= (Space2D.Vector) sim.space.makeVector();
			impulsefactor = 10;
			//impulsefactor = 50;
		}
    
		public void draw(java.awt.Graphics g, int[] origin, double scale)
		{
			if (!visible) return;
			g.setColor(Color.red);
			int x = origin[0]+(int)( targetpos.x(0) * scale*display.getToPixels());
			int y = origin[1]+(int)( targetpos.x(1) * scale*display.getToPixels());
			int xp = origin[0] + (int)(mousepos.x(0) *scale*display.getToPixels());
			int yp = origin[1] + (int)(mousepos.x(1) *scale*display.getToPixels());
			int radius = 25;
			g.drawOval(x-radius, y-radius, 2*radius, 2*radius);
			g.drawLine(x, y, xp, yp);
//			display.repaint();
		}
		public void setImpulseFactor(double d){
			impulsefactor =d;
		}
		public double getImpulseFactor(){
			return impulsefactor;
		}
		public void setAtom(Atom a){
			target = a;
		}
		public Atom getAtom(){
			return target;
		}
		public Component graphic(Object obj) {return new JPanel();}
    
    
		//  Start of DisplayPhaseListeners //
    
		private class Kicker implements DisplayPhaseListener {
			public void displayPhaseAction(DisplayPhaseEvent dpe) {
				java.awt.event.MouseEvent event =  dpe.getMouseEvent();
				targetpos = (Space2D.Vector) target.coord.position();
				impulse.E(dpe.point());
				impulse.ME(targetpos);
				impulse.TE(-impulsefactor*Math.exp(impulsefactor1*impulse.squared()));
				mousepos.E(dpe.point());
				if (event.getID() == java.awt.event.MouseEvent.MOUSE_PRESSED)
					visible = true;
				else if (event.getID() == java.awt.event.MouseEvent.MOUSE_RELEASED)
				{
					visible = false;}
					display.getPhase().integrator().pause();
					target.coord.accelerateBy(impulse);
					display.getPhase().integrator().unPause();
					display.getPhase().integrator().reset();
				
				display.repaint();
			}
		}
        
	}
    
}