package etomica.graphics;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.Panel;
import java.awt.TextField;

import etomica.atom.Atom;
import etomica.atom.AtomAgentManager;
import etomica.atom.AtomFilter;
import etomica.atom.AtomLeaf;
import etomica.atom.AtomTypeSphere;
import etomica.atom.AtomAgentManager.AgentIterator;
import etomica.atom.AtomAgentManager.AgentSource;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.integrator.IntegratorIntervalEvent;
import etomica.integrator.IntegratorIntervalListener;
import etomica.space.Vector;
import g3dsys.control.G3DSys;

//TODO: rewrite doPaint and drawAtom

public class DisplayPhaseCanvasG3DSys extends DisplayCanvas
	implements AgentSource, IntegratorIntervalListener {

	private TextField scaleText = new TextField();
	private final AtomIteratorLeafAtoms atomIterator = new AtomIteratorLeafAtoms();

	//will handle all actual drawing
	private G3DSys gsys;
	
	private AtomAgentManager aam;

	public DisplayPhaseCanvasG3DSys(DisplayPhase _phase) {
		//old stuff
		scaleText.setVisible(true);
		scaleText.setEditable(false);
		scaleText.setBounds(0,0,100,50);
		displayPhase = _phase;

		//init G3DSys
		Panel p = new Panel();
		this.setLayout(new java.awt.GridLayout());
		p.setLayout(new java.awt.GridLayout());
		p.setSize(800,800);
		this.add(p);
		gsys = new G3DSys(p);
		gsys.addFig(G3DSys.BOX, Color.GREEN, 0,0,0, 0);
		
		//init AtomAgentManager, to sync G3DSys and Etomica models
		//this automatically adds the atoms
		aam = new AtomAgentManager(this,displayPhase.getPhase());

		//need to init atoms here: get iterator, add figures, etc.
//		atomIterator.setPhase(displayPhase.getPhase());
//		atomIterator.reset();
//		while(atomIterator.hasNext()) {
//
//			AtomLeaf a = (AtomLeaf) atomIterator.nextAtom();
//			double[] coords = new double[3];
//			a.coord.position().assignTo(coords);
//			System.out.println("adding atom at "+coords[0]+","+coords[1]+","+coords[2]);
//			makeAgent(a);
//			gsys.addFig(G3DSys.Figures.BALL, org.jmol.g3d.Graphics3D.RED,
//			(float)coords[0], (float)coords[1], (float)coords[2], 1);
//		}
//		java.util.Random r = new java.util.Random();
//
//		for(int i=0; i<25000; i++) {
//			gsys.addFig(G3DSys.Figures.BALL, org.jmol.g3d.Graphics3D.BLUE,
//					20*r.nextFloat(), 20*r.nextFloat(), 20*r.nextFloat(), 1);
//		}
		gsys.refresh();
	}

	/**
	 * Sets the size of the display to a new value and scales the image so that
	 * the phase fits in the canvas in the same proportion as before.
	 */
	public void scaleSetSize(int width, int height) {
		if(getBounds().width * getBounds().height != 0) {  //reset scale based on larger size change
			double ratio1 = (double)width/(double)getBounds().width;
			double ratio2 = (double)height/(double)getBounds().height;
			double factor = Math.min(ratio1, ratio2);
			//        double factor = (Math.abs(Math.log(ratio1)) > Math.abs(Math.log(ratio2))) ? ratio1 : ratio2;
			displayPhase.setScale(displayPhase.getScale()*factor);
			setSize(width, height);
		}
	}

	//Override superclass methods for changing size so that scale is reset with any size change  
	// this setBounds is ultimately called by all other setSize, setBounds methods
	public void setBounds(int x, int y, int width, int height) {
		if(width == 0 || height == 0) return;
		super.setBounds(x,y,width,height);
		createOffScreen(width,height);
	}

	Vector vec2;  

	public void doPaint(Graphics g) {
        ColorScheme colorScheme = displayPhase.getColorScheme();
        AtomFilter atomFilter = displayPhase.getAtomFilter();
		atomIterator.setPhase(displayPhase.getPhase());
		atomIterator.reset();
		while(atomIterator.hasNext()) {
			AtomLeaf a = (AtomLeaf) atomIterator.nextAtom();
            if (!(a.type instanceof AtomTypeSphere)) continue;
			double[] coords = new double[3];
			a.coord.position().assignTo(coords);
			Long id = (Long)aam.getAgent(a);
            double diameter = ((AtomTypeSphere)a.type).getDiameter();
			gsys.modFig(id.longValue(), colorScheme.getAtomColor(a),
					(float)coords[0], (float)coords[1], (float)coords[2], (float)diameter,
                    atomFilter.accept(a));
		}
		gsys.fastRefresh();
		
	}

	/* ******************************************************
	 * AgentSource methods
	 * ******************************************************/
	public Class getAgentClass() {
		try {
			return java.lang.ClassLoader.getSystemClassLoader().loadClass("java.lang.Long");
		} catch(ClassNotFoundException cnfe) {
			System.out.println("OH NO! No Long type available.");
			return null;
		}
	}
	
	public Object makeAgent(Atom a) {
		if ( !(a instanceof AtomLeaf) || !(a.type instanceof AtomTypeSphere)) return null;
		double[] coords = new double[3];
		((AtomLeaf)a).coord.position().assignTo(coords);

        double diameter = ((AtomTypeSphere)a.type).getDiameter();
		long l = gsys.addFig(G3DSys.BALL, displayPhase.getColorScheme().getAtomColor((AtomLeaf)a),
				(float)coords[0], (float)coords[1], (float)coords[2], (float)diameter);
		//System.out.println("added atom "+l+" at "+coords[0]+","+coords[1]+","+coords[2]);
		//System.out.println("figs now: "+gsys.getFigs().length);
//		repaint();
//		gsys.fastRefresh();
		return new Long(l);
	}

	public void releaseAgent(Object agent, Atom atom) {
		//System.out.println("removed atom "+(Long)agent);
		gsys.removeFig(((Long) agent).longValue());
		//System.out.println("figs left: "+gsys.getFigs().length);
//		repaint();
//		gsys.fastRefresh();
	}

	/* *******************************************
	 * Integrator listener methods
	 * *******************************************/
	
	public int getPriority() {
		// TODO Auto-generated method stub
		return 0;
	}

	public void intervalAction(IntegratorIntervalEvent evt) {
		AgentIterator ai = aam.makeIterator();
		while(ai.hasNext()) {
			Long l = (Long) ai.next();
			System.out.println("got atom "+l);
		}
	}
	
}
