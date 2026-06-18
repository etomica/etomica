/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graphics;

import etomica.action.controller.Controller;
import etomica.atom.*;
import etomica.atom.AtomLeafAgentManager.AgentSource;
import etomica.box.Box;
import etomica.math.geometry.LineSegment;
import etomica.math.geometry.Plane;
import etomica.math.geometry.Polytope;
import etomica.space.Boundary;
import etomica.space.IOrientation;
import etomica.space.Vector;
import etomica.space3d.IOrientationFull3D;
import etomica.space3d.Space3D;
import etomica.units.Pixel;
import etomica.util.Arrays;
import g3dsys.control.G3DSys;
import g3dsys.images.*;
import org.jmol.util.Colix;
import org.jmol.util.Point3f;

import java.awt.*;
import java.util.HashMap;
import java.util.Map;

public class DisplayBoxCanvasG3DSys extends DisplayCanvas implements
        AgentSource<Ball> {


    protected final Map<AtomType, OrientedSite[]> atomTypeOrientedManager;
    private final double[] coords;
    protected AtomLeafAgentManager<Ball> aam;
    protected LineSegment[] lines;
    protected Line[] lineFigures;
    protected AtomLeafAgentManager<Ball[]> aamOriented;
    protected Vector rMin, rMax;
    // will handle all actual drawing
    protected G3DSys gsys;
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
    private Vector[] planeIntersections;
    private Vector work, work2, work3;
    private double[] planeAngles;
    private java.util.ArrayList<Object[]> pendingBonds = new java.util.ArrayList<Object[]>();

    public DisplayBoxCanvasG3DSys(DisplayBox _box, Controller controller) {
		super(controller);
		displayBox = _box;

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
        panel.setSize(2000, 1600);
        this.add(panel);
        coords = new double[3];
        gsys = new G3DSys(panel);
        setBackgroundColor(Color.BLACK);
        setBoundaryFrameColor(Color.WHITE);
        setPlaneColor(Color.YELLOW);
        // init AtomAgentManager, to sync G3DSys and Etomica models
        // this automatically adds the atoms
        aam = new AtomLeafAgentManager<Ball>(this, displayBox.getBox());
		OrientedAgentSource oas = new OrientedAgentSource();
		aamOriented = new AtomLeafAgentManager<Ball[]>(oas, displayBox.getBox());
		atomTypeOrientedManager = new HashMap<>();

        planes = new Plane[0];
        planeTriangles = new Triangle[0][0];
        planeIntersections = new Vector[0];
        lines = new LineSegment[0];
        lineFigures = new Line[0];
        planeAngles = new double[0];
        work = Space3D.getInstance().makeVector();
        work2 = Space3D.getInstance().makeVector();
        work3 = Space3D.getInstance().makeVector();

        pixel = new Pixel();
	}

	public G3DSys getG3DSys() {
	    return gsys;
	}

	/**
	 * Sets the display bounding box.  Atoms outside the box are not draw.  The
	 * boundary lines are snapped inside the bounding box if they are outside.
	 */
    public void setBoundingBox(Vector rMin, Vector rMax) {
        this.rMin = rMin;
        this.rMax = rMax;
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
     * Gets the background color of the display box canvas.
     *
     * @return Color : Current color of background
     */
    public Color getBackgroundColor() {
        return backgroundColor;
    }

	@Override
	public void setBackground(Color bg) {
		super.setBackground(bg);
		setBackgroundColor(bg);
	}

	/**
	 * Sets the background color of the display box canvas.
	 * @param color : color to set background to
	 */
	public void setBackgroundColor(Color color) {
		backgroundColor = color;
		if(gsys!=null) gsys.setBGColor(color);
		if(panel!=null) panel.setBackground(color);
    }

    /**
     * Gets the color of box boundary.
     * @return Color : Current color of box boundary
     */
    public Color getBoundaryFrameColor() {
        return boundaryFrameColor;
    }

	/**
	 * Sets the color of the box boundary.
	 * @param color : color to set box boundary
	 */
	public void setBoundaryFrameColor(Color color) {
		boundaryFrameColor = color;
		oldPolytope = null;
	}

	/**
	 * Gets the color of the plane.
	 * @return Color : Current color of plane
	 */
    public Color getPlaneColor() {
        return planeColor;
    }

    /**
     * Sets the color of the plane.
     *
     * @param color : color to set plane
     */
    public void setPlaneColor(Color color) {
        planeColor = color;
    }

	public void removeObjectByBox(Box p) {

		// Remove old box atoms
		IAtomList leafList = p.getLeafList();
		int nLeaf = leafList.size();
		for (int iLeaf = 0; iLeaf < nLeaf; iLeaf++) {
			IAtom a = leafList.get(iLeaf);
			if (a == null)
				continue;
			Ball ball = (Ball) aam.getAgent(a);
			if (ball == null) {
				continue;
			}
			gsys.removeFig(ball);
		}
	}

	/**
	 * refreshAtomAgentMgr() - sets the new atom manager based upon the box.
	 *    Would only need to be called if it's DisplayBoxs' box has changed.
	 *
	 */
	public void refreshAtomAgentMgr() {

		// Set new atom manager
        aam = new AtomLeafAgentManager<>(this, displayBox.getBox());
        aamOriented = new AtomLeafAgentManager<>(a -> null, displayBox.getBox());
		initialOrient = true;
	}

	/**
	 * helper method that returns the position r if there are no bounds or
	 * if r is inside the bounds or (if r is outside the bounds) returns the
	 * position in the bounds nearest r.
	 */
	protected double rBound(double r, int i) {
	    if (rMin == null) return r;
        return r < rMin.getX(i) ? rMin.getX(i) : (r > rMax.getX(i) ? rMax.getX(i) : r);
	}

	public void doPaint(Graphics g) {

		// handle pending bond addition requests
		if (pendingBonds.size() > 0) {
			for (int i = 0; i < pendingBonds.size(); i++) {
				Object[] o = pendingBonds.get(i);
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

		AtomTest atomTest = displayBox.getAtomTestDoDisplay();
		if (atomTest instanceof AtomTestCollective) {
			((AtomTestCollective) atomTest).resetTest();
        }
        ColorScheme colorScheme = displayBox.getColorScheme();
		if (colorScheme instanceof ColorSchemeCollective) {
			((ColorSchemeCollective) colorScheme).colorAllAtoms();
		}

		DiameterHash diameterHash = displayBox.getDiameterHash();

		IAtomList leafList = displayBox.getBox().getLeafList();
		int nLeaf = leafList.size();

		for (int iLeaf = 0; iLeaf < nLeaf; iLeaf++) {
		    IAtom a = null;
		    Ball ball = null;
		    try {
		        a = leafList.get(iLeaf);
	            if (a == null)
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
			boolean drawable = atomTest == null || atomTest.test(a);
            if (drawable && rMin != null) {
			    for (int i=0; i<rMin.getD(); i++) {
			        double x = a.getPosition().getX(i);
			        if (x < rMin.getX(i) || x > rMax.getX(i)) {
			            drawable = false;
			            break;
			        }
			    }
			}
			ball.setDrawable(drawable);
			if (drawable) {
				a.getPosition().assignTo(coords);
				float diameter = (float) diameterHash.getDiameter(a);
				// default diameter
				if (diameter == -1) diameter = 1;
				ball.setColor(G3DSys.getColix(colorScheme.getAtomColor(a)));
				ball.setD(diameter);
				ball.setX((float) coords[0]);
				ball.setY((float) coords[1]);
				ball.setZ((float) coords[2]);
			}

			OrientedSite[] sites = atomTypeOrientedManager.get(a.getType());
			if (sites != null) {
			    Ball[] ballSites = aamOriented.getAgent(a);
			    if (ballSites == null) {
					if (!drawable) continue;
					ballSites = new Ball[sites.length];
		            for (int j=0; j<sites.length; j++) {
		                ballSites[j] = new Ball(gsys, G3DSys.getColix(sites[j].color), 0, 0, 0, (float)sites[j].diameter);
		                gsys.addFig(ballSites[j]);
		            }
		            aamOriented.setAgent(a, ballSites);
			    } else {
					for (int i = 0; i < ballSites.length; i++) {
						ballSites[i].setDrawable(drawable);
					}
					if (!drawable) continue;
				}
				IOrientation orientation = ((IAtomOriented)a).getOrientation();
			    Vector direction1 = orientation.getDirection();
			    Vector direction2 = null;
			    if (orientation instanceof IOrientationFull3D) {
			        direction2 = ((IOrientationFull3D)orientation).getSecondaryDirection();
	                work2.E(direction1);
	                work2.XE(direction2);
			    }

			    for (int j=0; j<sites.length; j++) {
			        work.E(a.getPosition());
			        work.PEa1Tv1(sites[j].coord, direction1);
			        if (sites[j] instanceof OrientedFullSite) {
			            work.PEa1Tv1(((OrientedFullSite)sites[j]).coord2, direction2);
			            work.PEa1Tv1(((OrientedFullSite)sites[j]).coord3, work2);
			        }
			        work.assignTo(coords);
			        ballSites[j].setX((float) coords[0]);
			        ballSites[j].setY((float) coords[1]);
			        ballSites[j].setZ((float) coords[2]);
			    }
			}
		}

        for (int i=0; i<lines.length; i++) {
            updateLine(i);
        }

        for (int i=0; i<planes.length; i++) {
            updatePlane(i);
        }

		Boundary boundary = displayBox.getBox().getBoundary();

		// Do not draw bounding box around figure if the boundary
		// is not an etomica.space.Boundary

		Polytope polytope = boundary.getShape();
		if (polytope != oldPolytope) {

			// send iterator to g3dsys
			gsys.setBoundaryVectorsIterator(wrapIndexIterator(boundary
					.getIndexIterator()));

			if (polytopeLines != null) {
				for (int i = 0; i < polytopeLines.length; i++) {
					gsys.removeFig(polytopeLines[i]);
				}
			}
			LineSegment[] boundaryLines = polytope.getEdges();
			polytopeLines = new Line[boundaryLines.length];
			for (int i = 0; i < boundaryLines.length; i++) {
				Vector[] vertices = boundaryLines[i].getVertices();
				float v0x = (float)rBound(vertices[0].getX(0), 0);
				float v0y = (float)rBound(vertices[0].getX(1), 1);
				float v0z = (float)rBound(vertices[0].getX(2), 2);
				float v1x = (float)rBound(vertices[1].getX(0), 0);
				float v1y = (float)rBound(vertices[1].getX(1), 1);
				float v1z = (float)rBound(vertices[1].getX(2), 2);
				polytopeLines[i] = new Line(gsys, G3DSys.getColix(boundaryFrameColor),
					  Point3f.new3(v0x, v0y, v0z), Point3f.new3(v1x, v1y, v1z));
				if (displayBox.getShowBoundary() && drawBoundary > DRAW_BOUNDARY_NONE) {
					gsys.addFig(polytopeLines[i]);
				}
			}
			oldPolytope = polytope;
		} else {
			LineSegment[] boundaryLines = polytope.getEdges();
			for (int i = 0; i < boundaryLines.length; i++) {
				Vector[] vertices = boundaryLines[i].getVertices();

				float v0x = (float)rBound(vertices[0].getX(0), 0);
				float v0y = (float)rBound(vertices[0].getX(1), 1);
				float v0z = (float)rBound(vertices[0].getX(2), 2);
				float v1x = (float)rBound(vertices[1].getX(0), 0);
				float v1y = (float)rBound(vertices[1].getX(1), 1);
				float v1z = (float)rBound(vertices[1].getX(2), 2);
				polytopeLines[i].setStart(v0x, v0y, v0z);
				polytopeLines[i].setEnd(v1x, v1y, v1z);

				if ((!displayBox.getShowBoundary() || drawBoundary == DRAW_BOUNDARY_NONE) && boundaryDisplayed) {
					gsys.removeFig(polytopeLines[i]);
				} else if ((displayBox.getShowBoundary() && drawBoundary > DRAW_BOUNDARY_NONE) && !boundaryDisplayed) {
					gsys.addFig(polytopeLines[i]);
				}
			}
		}

		boundaryDisplayed = displayBox.getShowBoundary() && drawBoundary > DRAW_BOUNDARY_NONE;


		// set boundary vectors for image shell
		int n=0;
		for (int i=0; i<3; i++) {
			if (boundary.getPeriodicity(i)) {
				n++;
			}
		}
		double[] dvecs = new double[n * 3]; // assuming
														// 3-dimensional vectors
		int j = 0;
		for (int i = 0; i < 3; i++) {
			Vector v = boundary.getEdgeVector(i);
			if (!boundary.getPeriodicity(i)) {
				continue;
			}
			dvecs[j * 3] = v.getX(0);
			dvecs[j * 3 + 1] = v.getX(1);
			dvecs[j * 3 + 2] = v.getX(2);
			j++;
		}
		gsys.setBoundaryVectors(dvecs);

		Vector bounds = boundary.getBoxSize();
		gsys.setBoundingBox((float) (-bounds.getX(0) * 0.5),
				(float) (-bounds.getX(1) * 0.5), (float) (-bounds.getX(2) * 0.5),
				(float) (bounds.getX(0) * 0.5), (float) (bounds.getX(1) * 0.5),
				(float) (bounds.getX(2) * 0.5));

		// If displaying a new box, make it fit on the screen
		if(initialOrient == true) {
			gsys.scaleFitToScreen();
			initialOrient = false;
		}

		gsys.fastRefresh();
	}

    public void addLine(LineSegment newLine) {
        lines = (LineSegment[])Arrays.addObject(lines, newLine);
        Vector[] endpoints = newLine.getVertices();
        Point3f s = new Point3f();
        s.x = (float)endpoints[0].getX(0);
        s.y = (float)endpoints[0].getX(1);
        s.z = (float)endpoints[0].getX(2);
        Point3f e = new Point3f();
        e.x = (float)endpoints[1].getX(0);
        e.y= (float)endpoints[1].getX(1);
        e.z = (float)endpoints[1].getX(2);
        lineFigures = (Line[])Arrays.addObject(lineFigures, new Line(gsys, G3DSys.getColix(Color.WHITE), s, e));
        gsys.addFig(lineFigures[lineFigures.length-1]);
    }

    public void removeLine(LineSegment oldLine) {
        for (int i=0; i<lines.length; i++) {
            if (lines[i] == oldLine) {
                gsys.removeFig(lineFigures[i]);
                lineFigures = (Line[])Arrays.removeObject(lineFigures, lineFigures[i]);
                lines = (LineSegment[])Arrays.removeObject(lines, oldLine);
                return;
            }
        }
        throw new RuntimeException("I don't know about that line");
    }

    protected synchronized void updateLine(int iLine) {
        Vector[] endpoints = lines[iLine].getVertices();
        Line myLine = lineFigures[iLine];
        Point3f s = myLine.getStart();
        s.x = (float)endpoints[0].getX(0);
        s.y = (float)endpoints[0].getX(1);
        s.z = (float)endpoints[0].getX(2);
        Point3f e = myLine.getEnd();
        e.x = (float)endpoints[1].getX(0);
        e.y = (float)endpoints[1].getX(1);
        e.z = (float)endpoints[1].getX(2);
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
    
    protected synchronized void updatePlane(int iPlane) {
        Boundary boundary = displayBox.getBox().getBoundary();
    	if(!(boundary instanceof Boundary)) {
    		throw new RuntimeException("Unable to drawPlane for a Boundary not a subclass of etomica.space.Boundary");
    	}
        Plane plane = planes[iPlane];
        Polytope polytope = ((Boundary)boundary).getShape();
        LineSegment[] boundaryLines = polytope.getEdges();
        int intersectionCount = 0;
        for (int i = 0; i < boundaryLines.length; i++) {
            Vector[] vertices = boundaryLines[i].getVertices();
            work.Ev1Mv2(vertices[1], vertices[0]);
            // this happens to do what we want
            double alpha = -plane.distanceTo(vertices[0]) /
                            (plane.distanceTo(work) - plane.getD());
            if (alpha >= 0 && alpha <= 1) {
                Vector newIntersection;
                if (planeIntersections.length == intersectionCount) {
                    newIntersection = Space3D.getInstance().makeVector();
                    planeIntersections = (Vector[])Arrays.addObject(planeIntersections, newIntersection);
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
                    Vector intersection = planeIntersections[i];
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
                    gsys, Colix.getColixTranslucent3(G3DSys.getColix(planeColor), true, 0.5f), new Point3f(), new Point3f(), new Point3f()));
            gsys.addFig(planeTriangles[iPlane][planeTriangles[iPlane].length-1]);
        }

        for (int i=0; i<intersectionCount; i++) {
            planeIntersections[i].PE(work);
        }

        for (int i=0; i<planeTriangles[iPlane].length; i++) {
            Triangle triangle = planeTriangles[iPlane][i];
            Point3f p = triangle.getVertex1();
            p.x = (float)planeIntersections[0].getX(0);
            p.y = (float)planeIntersections[0].getX(1);
            p.z = (float)planeIntersections[0].getX(2);
            p = triangle.getVertex2();
            p.x = (float)planeIntersections[i+1].getX(0);
            p.y = (float)planeIntersections[i+1].getX(1);
            p.z = (float)planeIntersections[i+1].getX(2);
            p = triangle.getVertex3();
            p.x = (float)planeIntersections[i+2].getX(0);
            p.y = (float)planeIntersections[i+2].getX(1);
            p.z = (float)planeIntersections[i+2].getX(2);
        }
    }

    /**
	 * Add a bond to the graphical display between the given pairs. The given
	 * bondType is used to decide how the bond should be drawn.
	 */
	public Object makeBond(IAtomList pair, Object bondType) {
		/*
		 * Ball objects here could be null if the bond is created before the
		 * atoms have been added. Check for this and store atoms locally in a
		 * list. In doPaint check list for pending additions and add them.
		 */
		// bondType is a potential right now
		// best to ignore it for now; all bonds are equal
		Ball ball0 = (Ball) aam.getAgent(pair.get(0));
		Ball ball1 = (Ball) aam.getAgent(pair.get(1));
		if (ball0 == null || ball1 == null) {
			System.out.println("NULL!!!");
			pendingBonds.add(new Object[] { ball0, ball1, bondType });
			return null;
		}

		// make a bond object (Figure)
		Figure f = null;
		if (bondType == null || !(bondType instanceof Color)) {
		    f = new Bond(gsys, ball0, ball1);
		}
		else {
		    short c = G3DSys.getColix((Color)bondType);
		    f = new Bond(gsys, ball0, ball1, c);
		}
		gsys.addFig(f);
        return f;
    }

    public void setOrientationSites(AtomTypeOriented atomType, OrientedSite[] sites) {
        atomTypeOrientedManager.put(atomType, sites);
		IAtomList leafList = displayBox.getBox().getLeafList();
		int nLeaf = leafList.size();

		for (int iLeaf = 0; iLeaf < nLeaf; iLeaf++) {
			IAtom a = leafList.get(iLeaf);
			if (a.getType() != atomType) continue;
			Ball[] balls = aamOriented.getAgent(a);
			if (balls == null) continue;
			for (int j = 0; j < balls.length; j++) {
				gsys.removeFig(balls[j]);
			}
			aamOriented.setAgent(a, null);
		}
	}
	
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
    public Ball makeAgent(IAtom a, Box agentBox) {
		a.getPosition().assignTo(coords);

		float diameter = (float) displayBox.getDiameterHash().getDiameter(a);
		if (diameter == -1) diameter = 1;
		Ball newBall = new Ball(gsys, G3DSys.getColix((displayBox
				.getColorScheme().getAtomColor(a))), (float) coords[0],
				(float) coords[1], (float) coords[2], diameter);
		gsys.addFig(newBall);

		AtomTest atomFilter = displayBox.getAtomTestDoDisplay();
		if (atomFilter instanceof AtomTestCollective) {
			((AtomTestCollective) atomFilter).resetTest();
        }
        boolean drawable = atomFilter == null || atomFilter.test(a);
        if (drawable && rMin != null) {
            for (int i=0; i<rMin.getD(); i++) {
                double x = a.getPosition().getX(i);
                if (x < rMin.getX(i) || x > rMax.getX(i)) {
                    drawable = false;
                    break;
                }
            }
        }
        newBall.setDrawable(drawable);

		return newBall;
	}

    public void releaseAgent(Ball agent, IAtom atom, Box agentBox) {
        gsys.removeFig(agent);
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
    public double getDepth() {
        return gsys.getDepthPercent();
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

    public static class OrientedSite {
        public final double coord;
        public final Color color;
        public final double diameter;

        public OrientedSite(double coord, Color color, double diameter) {
            this.coord = coord;
            this.color = color;
            this.diameter = diameter;
        }
    }

    public static class OrientedFullSite extends OrientedSite {
        public final double coord2, coord3;

        public OrientedFullSite(Vector coord, Color color, double diameter) {
            super(coord.getX(0), color, diameter);
            coord2 = coord.getX(1);
            coord3 = coord.getX(2);
        }
    }

    public class OrientedAgentSource implements AgentSource<Ball[]> {

        public Ball[] makeAgent(IAtom a, Box agentBox) {
            return null;
        }

        public void releaseAgent(Ball[] agent, IAtom atom, Box agentBox) {
            for (int i = 0; i < agent.length; i++) {
                gsys.removeFig(agent[i]);
            }
        }
    }

}
