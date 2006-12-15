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
import etomica.atom.AtomAgentManager.AgentSource;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import g3dsys.control.G3DSys;
import g3dsys.images.Ball;
import g3dsys.images.Box;
import g3dsys.images.Figure;

//TODO: rewrite doPaint and drawAtom

public class DisplayPhaseCanvasG3DSys extends DisplayCanvas
	implements AgentSource {

	private TextField scaleText = new TextField();
	private final AtomIteratorLeafAtoms atomIterator = new AtomIteratorLeafAtoms();

	//will handle all actual drawing
	private G3DSys gsys;
    private final double[] coords;
	
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
        coords = new double[3];
		gsys = new G3DSys(p);
		gsys.addFig(new Box(gsys, G3DSys.getColix(Color.GREEN)));
		
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

	public void doPaint(Graphics g) {
        ColorScheme colorScheme = displayPhase.getColorScheme();
        AtomFilter atomFilter = displayPhase.getAtomFilter();
		atomIterator.setPhase(displayPhase.getPhase());
		atomIterator.reset();
		while(atomIterator.hasNext()) {
			AtomLeaf a = (AtomLeaf) atomIterator.nextAtom();
            if (!(a.getType() instanceof AtomTypeSphere)) continue;
            Ball ball = (Ball)aam.getAgent(a);
            if (ball == null) {
                continue;
            }
            boolean drawable = atomFilter.accept(a);
            ball.setDrawable(drawable);
            if (!drawable) {
                continue;
            }
            a.getCoord().position().assignTo(coords);
            float diameter = (float)((AtomTypeSphere)a.getType()).getDiameter();
            ball.setColor(G3DSys.getColix(colorScheme.getAtomColor(a)));
            ball.setD(diameter);
            ball.setX((float)coords[0]);
            ball.setY((float)coords[1]);
            ball.setZ((float)coords[2]);
		}
		gsys.fastRefresh();
		
	}

	/* ******************************************************
	 * AgentSource methods
	 * ******************************************************/
	public Class getAgentClass() {
	    return Figure.class;
	}
	
	public Object makeAgent(Atom a) {
		if ( !(a instanceof AtomLeaf) || !(a.getType() instanceof AtomTypeSphere)) return null;
		((AtomLeaf)a).getCoord().position().assignTo(coords);

        float diameter = (float)((AtomTypeSphere)a.getType()).getDiameter();
        Ball newBall = new Ball(gsys, G3DSys.getColix((displayPhase.getColorScheme().getAtomColor((AtomLeaf)a))),
                (float)coords[0], (float)coords[1], (float)coords[2], diameter);
        gsys.addFig(newBall);
        return newBall;
	}

	public void releaseAgent(Object agent, Atom atom) {
		gsys.removeFig((Figure) agent);
	}
}
