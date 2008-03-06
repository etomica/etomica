package etomica.graphics;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.Panel;
import java.awt.TextField;

import javax.vecmath.Point3f;

import org.jmol.g3d.Graphics3D;

import etomica.api.IAtom;
import etomica.api.IAtomPositioned;
import etomica.api.IAtomSet;
import etomica.api.IBoundary;
import etomica.api.IBox;
import etomica.api.IVector;
import etomica.atom.AtomFilter;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomTypeSphere;
import etomica.atom.AtomAgentManager.AgentSource;
import etomica.space.Boundary;
import etomica.math.geometry.LineSegment;
import etomica.math.geometry.Plane;
import etomica.math.geometry.Polytope;
import etomica.space.Space;
import etomica.space3d.IVector3D;
import etomica.util.Arrays;
import g3dsys.control.G3DSys;
import g3dsys.images.Ball;
import g3dsys.images.Bond;
import g3dsys.images.Figure;
import g3dsys.images.Line;
import g3dsys.images.Triangle;

public class DisplayBoxCanvasG3DSys extends DisplayCanvas implements
		AgentSource, BondManager {

	private TextField scaleText = new TextField();

	// will handle all actual drawing
	private G3DSys gsys;
	private final double[] coords;

	private AtomLeafAgentManager aam;

	private Polytope oldPolytope;
	private Line[] polytopeLines;
	private boolean boundaryDisplayed = false;
	private Color backgroundColor;
	private Color boundaryFrameColor;
	private Color planeColor;
	private Panel panel = null;
	private boolean initialOrient = false;
    private Plane[] planes;
    private Triangle[][] planeTriangles;
    private IVector3D[] planeIntersections;
    private IVector3D work, work2, work3;
    private double[] planeAngles;
    private final Space space;

	public DisplayBoxCanvasG3DSys(DisplayBox _box, Space _space) {
		// old stuff
		scaleText.setVisible(true);
		scaleText.setEditable(false);
		scaleText.setBounds(0, 0, 100, 50);
		displayBox = _box;
		space = _space;

		// init G3DSys
		// adding JPanel flickers, Panel does not. Nobody knows why.
		/*
		 * Set visible false here to be toggled later; seems to fix the
		 * 'sometimes gray' bug
		 */
		// this.setVisible(false); // to be set visible later by
		// SimulationGraphic
		panel = new Panel();
		this.setLayout(new java.awt.GridLayout());
		panel.setLayout(new java.awt.GridLayout());
		panel.setSize(800, 800);
		this.add(panel);
		coords = new double[3];
		gsys = new G3DSys(panel);
		setBackgroundColor(Color.BLACK);
		setBoundaryFrameColor(Color.WHITE);
		setPlaneColor(Color.YELLOW);
		// init AtomAgentManager, to sync G3DSys and Etomica models
		// this automatically adds the atoms
		aam = new AtomLeafAgentManager(this, displayBox.getBox(), false);

		planes = new Plane[0];
        planeTriangles = new Triangle[0][0];
        planeIntersections = new IVector3D[0];
        planeAngles = new double[0];
        work = (IVector3D)space.makeVector();
        work2 = (IVector3D)space.makeVector();
        work3 = (IVector3D)space.makeVector();
	}

	/**
	 * Sets the size of the display to a new value and scales the image so that
	 * the box fits in the canvas in the same proportion as before.
	 */
	public void scaleSetSize(int width, int height) {
		if (getBounds().width * getBounds().height != 0) { // reset scale based
															// on larger size
															// change
			double ratio1 = (double) width / (double) getBounds().width;
			double ratio2 = (double) height / (double) getBounds().height;
			double factor = Math.min(ratio1, ratio2);
			// double factor = (Math.abs(Math.log(ratio1)) >
			// Math.abs(Math.log(ratio2))) ? ratio1 : ratio2;
			displayBox.setScale(displayBox.getScale() * factor);
			setSize(width, height);
		}
	}

	// Override superclass methods for changing size so that scale is reset with
	// any size change
	// this setBounds is ultimately called by all other setSize, setBounds
	// methods
	public void setBounds(int x, int y, int width, int height) {
		if (width <= 0 || height <= 0)
			return;
		super.setBounds(x, y, width, height);
		createOffScreen(width, height);
	}

	/**
	 * Sets the background color of the display box canvas.
	 * @param Color : color to set background to
	 */
	public void setBackgroundColor(Color color) {
		backgroundColor = color;
		gsys.setBGColor(color);
		panel.setBackground(color);
	}

	/**
	 * Gets the background color of the display box canvas.
	 * @return Color : Current color of background
	 */
	public Color getBackgroundColor() {
		return backgroundColor;
	}

	/**
	 * Sets the color of the box boundary.
	 * @param Color : color to set box boundary
	 */
	public void setBoundaryFrameColor(Color color) {
		boundaryFrameColor = color;
		oldPolytope = null;
	}

	/**
	 * Gets the color of box boundary.
	 * @return Color : Current color of box boundary
	 */
	public Color getBoundaryFrameColor() {
		return boundaryFrameColor;
	}

	/**
	 * Sets the color of the plane.
	 * @param Color : color to set plane
	 */
	public void setPlaneColor(Color color) {
		planeColor = color;
	}

	/**
	 * Gets the color of the plane.
	 * @return Color : Current color of plane
	 */
	public Color getPlaneColor(Color color) {
		return planeColor;
	}

	public void removeObjectByBox(IBox p) {

		// Remove old box atoms
		IAtomSet leafList = p.getLeafList();
		int nLeaf = leafList.getAtomCount();
		for (int iLeaf = 0; iLeaf < nLeaf; iLeaf++) {
			IAtomPositioned a = (IAtomPositioned) leafList.getAtom(iLeaf);
			if (a == null || !(a.getType() instanceof AtomTypeSphere))
				continue;
			Ball ball = (Ball) aam.getAgent(a);
			if (ball == null) {
				continue;
			} else {
				gsys.removeFig(ball);
			}
		}
	}

	/**
	 * refreshAtomAgentMgr() - sets the new atom manager based upon the box.
	 *    Would only need to be called if it's DisplayBoxs' box has changed.
	 *
	 */
	public void refreshAtomAgentMgr() {

		// Set new atom manager
		aam = null;
		aam = new AtomLeafAgentManager(this, displayBox.getBox(), false);
		initialOrient = true;
	}

	public void doPaint(Graphics g) {

		// handle pending bond addition requests
		if (pendingBonds.size() > 0) {
			for (int i = 0; i < pendingBonds.size(); i++) {
				Object[] o = (Object[]) pendingBonds.get(i);
				Ball ball0 = (Ball) o[0];
				Ball ball1 = (Ball) o[1];
				// can't do anything with bondType for now
				Figure f = new Bond(gsys, ball0, ball1);
				gsys.addFig(f);
			}
		}

/*
 UNCOMMENT THIS CODE IF APPLICATIONS START ADDING THEIR OWN DRAWABLES

        //do drawing of all drawing objects that have been added to the display
        for(Iterator iter=displayBox.getDrawables().iterator(); iter.hasNext(); ) {
            Drawable obj = (Drawable)iter.next();
            obj.draw(g, displayBox.getOrigin(), displayBox.getToPixels());
        }
*/

		ColorScheme colorScheme = displayBox.getColorScheme();
		AtomFilter atomFilter = displayBox.getAtomFilter();
		if (colorScheme instanceof ColorSchemeCollective) {
			((ColorSchemeCollective) colorScheme).colorAllAtoms();
		}

		IAtomSet leafList = displayBox.getBox().getLeafList();
		int nLeaf = leafList.getAtomCount();

		for (int iLeaf = 0; iLeaf < nLeaf; iLeaf++) {
		    IAtomPositioned a = null;
		    Ball ball = null;
		    try {
		        a = (IAtomPositioned) leafList.getAtom(iLeaf);
	            if (a == null || !(a.getType() instanceof AtomTypeSphere))
	                continue;
	            ball = (Ball) aam.getAgent(a);
		    }
		    catch (ArrayIndexOutOfBoundsException e) {
		        System.out.println("oops, array index out of bounds");
		        //atoms might have been removed on another thread
		        break;
		    }
            catch (IndexOutOfBoundsException e) {
                System.out.println("oops, index out of bounds");
                //atoms might have been removed on another thread
                break;
            }
			if (ball == null) {
				continue;
			}
			/*
			 * Atomfilter changes the drawable flag in spheres; bonds respect
			 * this and will not draw themselves either. Wireframe mode, on the
			 * other hand, tells G3DSys to ignore spheres entirely regardless of
			 * drawable flag. This makes it possible to filter bonds in
			 * wireframe mode as well.
			 */
			boolean drawable = atomFilter.accept(a);
			ball.setDrawable(drawable);
			if (!drawable) {
				continue;
			}
			a.getPosition().assignTo(coords);
			float diameter = (float) ((AtomTypeSphere) a.getType())
					.getDiameter();
			ball.setColor(G3DSys.getColix(colorScheme.getAtomColor(a)));
			ball.setD(diameter);
			ball.setX((float) coords[0]);
			ball.setY((float) coords[1]);
			ball.setZ((float) coords[2]);
		}

        for (int i=0; i<planes.length; i++) {
            drawPlane(i);
        }
        
		IBoundary boundary = displayBox.getBox().getBoundary();

		// Do not draw bounding box around figure if the boundary
		// is not an etomica.space.Boundary
		if(boundary instanceof Boundary) {

			Polytope polytope = ((Boundary)boundary).getShape();
			if (polytope != oldPolytope) {
	
				// force trunc. oct. to make vecs else null pointer exception
				boundary.getPeriodicVectors();
				// send iterator to g3dsys
				gsys.setBoundaryVectorsIterator(wrapIndexIterator((boundary
						.getIndexIterator())));
	
				if (polytopeLines != null) {
					for (int i = 0; i < polytopeLines.length; i++) {
						gsys.removeFig(polytopeLines[i]);
					}
				}
				LineSegment[] lines = polytope.getEdges();
				polytopeLines = new Line[lines.length];
				for (int i = 0; i < lines.length; i++) {
					IVector[] vertices = lines[i].getVertices();
					polytopeLines[i] = new Line(gsys, G3DSys
							.getColix(boundaryFrameColor), new Point3f(
							(float) vertices[0].x(0), (float) vertices[0].x(1),
							(float) vertices[0].x(2)), new Point3f(
							(float) vertices[1].x(0), (float) vertices[1].x(1),
							(float) vertices[1].x(2)));
					if (displayBox.getShowBoundary() == true) {
						gsys.addFig(polytopeLines[i]);
					}
				}
				oldPolytope = polytope;
			} else {
				LineSegment[] lines = polytope.getEdges();
				for (int i = 0; i < lines.length; i++) {
					IVector[] vertices = lines[i].getVertices();
					polytopeLines[i].setStart((float) vertices[0].x(0),
							(float) vertices[0].x(1), (float) vertices[0].x(2));
					polytopeLines[i].setEnd((float) vertices[1].x(0),
							(float) vertices[1].x(1), (float) vertices[1].x(2));
	
					if (displayBox.getShowBoundary() == false
							&& boundaryDisplayed == true) {
						gsys.removeFig(polytopeLines[i]);
					} else if (displayBox.getShowBoundary() == true
							&& boundaryDisplayed == false) {
						gsys.addFig(polytopeLines[i]);
					}
				}
			}
	
			if (displayBox.getShowBoundary() == false) {
				boundaryDisplayed = false;
			} else {
				boundaryDisplayed = true;
			}

		}

		// set boundary vectors for image shell
		IVector[] vecs = boundary.getPeriodicVectors();
		double[] dvecs = new double[vecs.length * 3]; // assuming
														// 3-dimensional vectors
		for (int i = 0; i < vecs.length; i++) {
			if (vecs[i] == null)
				continue;
			dvecs[i * 3] = vecs[i].x(0);
			dvecs[i * 3 + 1] = vecs[i].x(1);
			dvecs[i * 3 + 2] = vecs[i].x(2);
		}
		gsys.setBoundaryVectors(dvecs);

		IVector bounds = boundary.getBoundingBox();
		gsys.setBoundingBox((float) (-bounds.x(0) * 0.5),
				(float) (-bounds.x(1) * 0.5), (float) (-bounds.x(2) * 0.5),
				(float) (bounds.x(0) * 0.5), (float) (bounds.x(1) * 0.5),
				(float) (bounds.x(2) * 0.5));

		// If displaying a new box, make it fit on the screen
		if(initialOrient == true) {
			gsys.scaleFitToScreen();
			initialOrient = false;
		}

		gsys.fastRefresh();
	}

    public void addPlane(Plane newPlane) {
        planes = (Plane[])Arrays.addObject(planes, newPlane);
        planeTriangles = (Triangle[][])Arrays.addObject(planeTriangles, new Triangle[0]);
    }
    
    public void removePlane(Plane oldPlane) {
        for (int i=0; i<planes.length; i++) {
            if (planes[i] == oldPlane) {
                for (int j=0; j<planeTriangles[i].length; j++) {
                    gsys.removeFig(planeTriangles[i][j]);
                }
                planeTriangles = (Triangle[][])Arrays.removeObject(planeTriangles, planeTriangles[i]);
                planes = (Plane[])Arrays.removeObject(planes, oldPlane);
                return;
            }
        }
        throw new RuntimeException("I don't know about that plane");
    }
        
    public synchronized void drawPlane(int iPlane) {
        IBoundary boundary = displayBox.getBox().getBoundary();
    	if(!(boundary instanceof Boundary)) {
    		throw new RuntimeException("Unable to drawPlane for a Boundary not a subclass of etomica.space.Boundary");
    	}
        Plane plane = planes[iPlane];
        Polytope polytope = ((Boundary)boundary).getShape();
        LineSegment[] lines = polytope.getEdges();
        int intersectionCount = 0;
        for (int i = 0; i < lines.length; i++) {
            IVector[] vertices = lines[i].getVertices();
            work.Ev1Mv2(vertices[1], vertices[0]);
            // this happens to do what we want
            double alpha = -plane.distanceTo(vertices[0]) / 
                            (plane.distanceTo(work) - plane.getD());
            if (alpha >= 0 && alpha <= 1) {
                IVector newIntersection;
                if (planeIntersections.length == intersectionCount) {
                    newIntersection = space.makeVector();
                    planeIntersections = (IVector3D[])Arrays.addObject(planeIntersections, newIntersection);
                }
                else {
                    newIntersection = planeIntersections[intersectionCount];
                }
                intersectionCount++;
                newIntersection.E(vertices[0]);
                newIntersection.PEa1Tv1(alpha, work);
            }
        }
        if (intersectionCount < 3) {
            for (int i=0; i<planeTriangles[iPlane].length; i++) {
                gsys.removeFig(planeTriangles[iPlane][i]);
            }
            planeTriangles[iPlane] = new Triangle[0];
            return;
        }
        
        //find the center of the polygon
        work.E(0);
        for (int i=0; i<intersectionCount; i++) {
            work.PE(planeIntersections[i]);
        }
        work.TE(1.0/intersectionCount);
        
        // convert the vertices to be vectors from the center
        // we'll switch back later
        for (int i=0; i<intersectionCount; i++) {
            planeIntersections[i].ME(work);
        }

        if (planeAngles.length < intersectionCount-1) {
            planeAngles = new double[intersectionCount-1];
        }
        work2.E(planeIntersections[0]);
        work2.XE(planeIntersections[1]);
        for (int i=1; i<intersectionCount; i++) {
            // If you understood this without reading this comment, I'll be
            // impressed.  The purpose here is to put the array of
            // intersections in order such that they form a polygon and lines
            // drawn between consecutive points don't intersect.  So we
            // calculate the angle between the first polygon point (which is
            // arbitrary), the center and each other point.  And we sort the
            // points by that angle.  We check the cross product so we can
            // distinguish 30 degrees from 330 degrees.
            double dot = planeIntersections[0].dot(planeIntersections[i]);
            double angle = Math.acos(dot / Math.sqrt(planeIntersections[0].squared() * planeIntersections[i].squared()));
            work3.E(planeIntersections[0]);
            work3.XE(planeIntersections[i]);
            // work2 dot work3 should be |work2|^2 or -|work2|^2.  Positive
            // indicates the angle is <180, negative indicates >180.
            if (work3.dot(work2) < 0) {
                angle = 2*Math.PI - angle;
            }
            boolean success = false;
            for (int j=1; j<i; j++) {
                if (angle < planeAngles[j-1]) {
                    // insert the i point at position j, shift existing points
                    IVector3D intersection = planeIntersections[i];
                    for (int k=i; k>j; k--) {
                        planeAngles[k-1] = planeAngles[k-2];
                        planeIntersections[k] = planeIntersections[k-1];
                    }
                    planeIntersections[j] = intersection;
                    planeAngles[j-1] = angle;
                    success = true;
                    break;
                }
            }
            if (!success) {
                planeAngles[i-1] = angle;
            }
        }

        // we need N-2 triangles
        while (intersectionCount < planeTriangles[iPlane].length+2) {
            Triangle triangle = planeTriangles[iPlane][planeTriangles[iPlane].length-1];
            gsys.removeFig(planeTriangles[iPlane][planeTriangles[iPlane].length-1]);
            planeTriangles[iPlane] = (Triangle[])Arrays.removeObject(planeTriangles[iPlane], triangle);
        }
        while (intersectionCount > planeTriangles[iPlane].length+2) {
            planeTriangles[iPlane] = (Triangle[])Arrays.addObject(planeTriangles[iPlane], new Triangle(
                    gsys, Graphics3D.getColixTranslucent(G3DSys.getColix(planeColor), true), new Point3f(), new Point3f(), new Point3f()));
            gsys.addFig(planeTriangles[iPlane][planeTriangles[iPlane].length-1]);
        }

        for (int i=0; i<intersectionCount; i++) {
            planeIntersections[i].PE(work);
        }

        for (int i=0; i<planeTriangles[iPlane].length; i++) {
            Triangle triangle = planeTriangles[iPlane][i];
            Point3f p = triangle.getVertex1();
            p.x = (float)planeIntersections[0].x(0);
            p.y = (float)planeIntersections[0].x(1);
            p.z = (float)planeIntersections[0].x(2);
            p = triangle.getVertex2();
            p.x = (float)planeIntersections[i+1].x(0);
            p.y = (float)planeIntersections[i+1].x(1);
            p.z = (float)planeIntersections[i+1].x(2);
            p = triangle.getVertex3();
            p.x = (float)planeIntersections[i+2].x(0);
            p.y = (float)planeIntersections[i+2].x(1);
            p.z = (float)planeIntersections[i+2].x(2);
        }
    }
    
	/**
	 * Add a bond to the graphical display between the given pairs. The given
	 * bondType is used to decide how the bond should be drawn.
	 */
	public Object makeBond(IAtomSet pair, Object bondType) {
		/*
		 * Ball objects here could be null if the bond is created before the
		 * atoms have been added. Check for this and store atoms locally in a
		 * list. In doPaint check list for pending additions and add them.
		 */
		// bondType is a potential right now
		// best to ignore it for now; all bonds are equal
		Ball ball0 = (Ball) aam.getAgent(pair.getAtom(0));
		Ball ball1 = (Ball) aam.getAgent(pair.getAtom(1));
		if (ball0 == null || ball1 == null) {
			System.out.println("NULL!!!");
			pendingBonds.add(new Object[] { ball0, ball1, bondType });
			return null;
		}

		// make a bond object (Figure)
		Figure f = new Bond(gsys, ball0, ball1);
		gsys.addFig(f);
		return f;
	}

	private java.util.ArrayList pendingBonds = new java.util.ArrayList();

	/**
	 * Removes the given bond from the graphical display. The bond must be an
	 * Object returned by the makeBond method.
	 */
	public void releaseBond(Object bond) {
		Figure figure = (Figure) bond;
		if (figure.getID() == -1) {
			throw new RuntimeException(figure + " has already been removed");
		}
		gsys.removeFig(figure);
	}

	/***************************************************************************
	 * AgentSource methods
	 **************************************************************************/
	public Class getAgentClass() {
		return Figure.class;
	}

	public Object makeAgent(IAtom a) {
		if (!(a.getType() instanceof AtomTypeSphere))
			return null;
		((IAtomPositioned) a).getPosition().assignTo(coords);

		float diameter = (float) ((AtomTypeSphere) a.getType()).getDiameter();
		Ball newBall = new Ball(gsys, G3DSys.getColix((displayBox
				.getColorScheme().getAtomColor(a))), (float) coords[0],
				(float) coords[1], (float) coords[2], diameter);
		gsys.addFig(newBall);
		return newBall;
	}

	public void releaseAgent(Object agent, IAtom atom) {
		gsys.removeFig((Figure) agent);
	}

	/**
	 * Set slab percentage
	 * 
	 * @param slab
	 *            the slab percentage to set
	 */
	public void setSlab(double slab) {
		gsys.setSlabPercent((int) slab);
	}

	/**
	 * Get depth percentage
	 * 
	 * @return returns current depth percentage
	 */
	public double getSlab() {
		return gsys.getSlabPercent();
	}

	/**
	 * Set depth percentage
	 * 
	 * @param depth
	 *            the depth percentage to set
	 */
	public void setDepth(double depth) {
		gsys.setDepthPercent((int) depth);
	}

	/**
	 * Get depth percentage
	 * 
	 * @return returns current depth percentage
	 */
	public double getDepth() {
		return gsys.getDepthPercent();
	}

	public void stopRotate() {
		gsys.stopRotation();
	}

	/**
	 * Wraps an etomica index iterator in an equivalent g3dsys interface for
	 * transport; removes g3dsys dependency from all but the etomica.graphics
	 * package.
	 * 
	 * @param iter
	 *            the etomica index iterator to wrap
	 * @return returns the g3dsys index iterator
	 */
	private g3dsys.control.IndexIterator wrapIndexIterator(
			etomica.lattice.IndexIteratorSizable iter) {

		final etomica.lattice.IndexIteratorSizable i = iter;

		return new g3dsys.control.IndexIterator() {

			private etomica.lattice.IndexIteratorSizable ii = i;

			public int getD() {
				return ii.getD();
			}

			public boolean hasNext() {
				return ii.hasNext();
			}

			public int[] next() {
				return ii.next();
			}

			public void reset() {
				ii.reset();
			}

			public void setSize(int[] size) {
				ii.setSize(size);
			}

			public boolean isLazySafe() {
				/*
				 * For now all boundaries are lazy-safe, including truncated
				 * octahedron. If this changes, check for that boundary type
				 * here (instanceof IndexIteratorSequentialFiltered, say, after
				 * making the class public) and use appropriate boolean.
				 */
				return true;
			}

		};

	}

}
