package etomica.graphics;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.Panel;
import java.awt.TextField;

import javax.vecmath.Point3f;

import etomica.atom.Atom;
import etomica.atom.AtomAgentManager;
import etomica.atom.AtomFilter;
import etomica.atom.AtomLeaf;
import etomica.atom.AtomTypeSphere;
import etomica.atom.AtomAgentManager.AgentSource;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.math.geometry.LineSegment;
import etomica.math.geometry.Polytope;
import etomica.space.Boundary;
import etomica.space.IVector;
import g3dsys.control.G3DSys;
import g3dsys.images.Ball;
import g3dsys.images.Figure;
import g3dsys.images.Line;

//TODO: rewrite doPaint and drawAtom

public class DisplayPhaseCanvasG3DSys extends DisplayCanvas
	implements AgentSource {

	private TextField scaleText = new TextField();
	private final AtomIteratorLeafAtoms atomIterator = new AtomIteratorLeafAtoms();

	//will handle all actual drawing
	private G3DSys gsys;
    private final double[] coords;
	
	private AtomAgentManager aam;
    
    private Polytope oldPolytope;
    private Line[] polytopeLines;

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
		
		//init AtomAgentManager, to sync G3DSys and Etomica models
		//this automatically adds the atoms
		aam = new AtomAgentManager(this, displayPhase.getPhase(), false);

//		gsys.refresh();
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
            a.getCoord().getPosition().assignTo(coords);
            float diameter = (float)((AtomTypeSphere)a.getType()).getDiameter();
            ball.setColor(G3DSys.getColix(colorScheme.getAtomColor(a)));
            ball.setD(diameter);
            ball.setX((float)coords[0]);
            ball.setY((float)coords[1]);
            ball.setZ((float)coords[2]);
		}
        
        Boundary boundary = displayPhase.getPhase().getBoundary();
        Polytope polytope = boundary.getShape();
        if (polytope != oldPolytope) {
            if (polytopeLines != null) {
                for (int i=0; i<polytopeLines.length; i++) {
                    gsys.removeFig(polytopeLines[i]);
                }
            }
            LineSegment[] lines = polytope.getEdges();
            polytopeLines = new Line[lines.length];
            for (int i=0; i<lines.length; i++) {
                IVector[] vertices = lines[i].getVertices();
                polytopeLines[i] = new Line(gsys, G3DSys.getColix(Color.WHITE), 
                        new Point3f((float)vertices[0].x(0), (float)vertices[0].x(1), (float)vertices[0].x(2)), 
                        new Point3f((float)vertices[1].x(0), (float)vertices[1].x(1), (float)vertices[1].x(2)));
                gsys.addFig(polytopeLines[i]);
            }
            oldPolytope = polytope;
        }
        else {
            LineSegment[] lines = polytope.getEdges();
            for (int i=0; i<lines.length; i++) {
                IVector[] vertices = lines[i].getVertices();
                polytopeLines[i].setStart((float)vertices[0].x(0), (float)vertices[0].x(1), (float)vertices[0].x(2));
                polytopeLines[i].setEnd((float)vertices[1].x(0), (float)vertices[1].x(1), (float)vertices[1].x(2));
            }
        }
        IVector bounds = boundary.getBoundingBox();
        gsys.setBoundingBox((float)(-bounds.x(0)*0.5), (float)(-bounds.x(1)*0.5), (float)(-bounds.x(2)*0.5),
                            (float)( bounds.x(0)*0.5), (float)( bounds.x(1)*0.5), (float)( bounds.x(2)*0.5));
        
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
		((AtomLeaf)a).getCoord().getPosition().assignTo(coords);

        float diameter = (float)((AtomTypeSphere)a.getType()).getDiameter();
        Ball newBall = new Ball(gsys, G3DSys.getColix((displayPhase.getColorScheme().getAtomColor((AtomLeaf)a))),
                (float)coords[0], (float)coords[1], (float)coords[2], diameter);
        gsys.addFig(newBall);
        return newBall;
	}

	public void releaseAgent(Object agent, Atom atom) {
		gsys.removeFig((Figure) agent);
	}
    
    
    public void setSlab(double slab) {
      gsys.setSlabPercent((int)slab);
    }
    public double getSlab() {
      return gsys.getSlabPercent();
    }
    public void setDepth(double depth) {
      gsys.setDepthPercent((int)depth);
    }
    public double getDepth() {
      return gsys.getDepthPercent();
    }
}
