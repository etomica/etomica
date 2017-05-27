/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graphics;

import etomica.space.Vector;
import etomica.math.geometry.LineSegment;
import etomica.math.geometry.Polytope;
import g3dsys.control.G3DSys;
import g3dsys.images.Line;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.Panel;

import org.jmol.util.Point3f;

public class DisplayPolytopeCanvasG3DSys extends DisplayCanvas {


	// will handle all actual drawing
	private G3DSys gsys;

	private Line[] polytopeLines;
	private Color backgroundColor;
	private Color boundaryFrameColor;
	private Panel panel = null;
	private boolean initialOrient = false;
    private DisplayPolytope displayPolytope;


	public DisplayPolytopeCanvasG3DSys(DisplayPolytope _box) {
	    super(null);
		displayPolytope = _box;

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
		panel.setSize(1600, 1600);
		this.add(panel);
		gsys = new G3DSys(panel);
		setBackgroundColor(Color.BLACK);
		setBoundaryFrameColor(Color.WHITE);
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
	 * @param color : color to set background to
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
	 * @param color : color to set box boundary
	 */
	public void setBoundaryFrameColor(Color color) {
		boundaryFrameColor = color;
		polytopeLines = null;
	}

	/**
	 * Gets the color of box boundary.
	 * @return Color : Current color of box boundary
	 */
	public Color getBoundaryFrameColor() {
		return boundaryFrameColor;
	}

	public void doPaint(Graphics g) {
	    Polytope polytope = displayPolytope.getPolytope();
			if (polytopeLines == null) {
	
			LineSegment[] lines = polytope.getEdges();
			polytopeLines = new Line[lines.length];
			for (int i = 0; i < lines.length; i++) {
				Vector[] vertices = lines[i].getVertices();
				polytopeLines[i] = new Line(gsys, G3DSys
						.getColix(boundaryFrameColor), Point3f.new3(
						(float) vertices[0].getX(0), (float) vertices[0].getX(1),
						(float) vertices[0].getX(2)), Point3f.new3(
						(float) vertices[1].getX(0), (float) vertices[1].getX(1),
						(float) vertices[1].getX(2)));
				gsys.addFig(polytopeLines[i]);
			}
		} else {
			LineSegment[] lines = polytope.getEdges();
			for (int i = 0; i < lines.length; i++) {
				Vector[] vertices = lines[i].getVertices();
				polytopeLines[i].setStart((float) vertices[0].getX(0),
						(float) vertices[0].getX(1), (float) vertices[0].getX(2));
				polytopeLines[i].setEnd((float) vertices[1].getX(0),
						(float) vertices[1].getX(1), (float) vertices[1].getX(2));
			}
		}

	    Vector bounds = displayPolytope.dimensions();
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
}
