package etomica.virial.cluster.graphics;

import java.awt.Canvas;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Image;
import java.awt.Point;
import java.util.Hashtable;

/**
 * 
 * Used from comphys.graphic.plot
 * A class which makes plotting a little more convenient than using
 * Java's Graphics methods.
 * 
 */
public abstract class Plot extends Canvas {

    protected int xPixels = 400;
    protected int yPixels = 400;

    /**
     * Constructor sets default size of the plot
     */
    public Plot () {
	setBackground(Color.white);
	setSize(xPixels, yPixels);
    }

    /**
     * Sets width and height of plot in pixels.
     *
     * @param x the plot width in pixels
     * @param y the plot height in pixels
     */
    public void setSize (int x, int y) {
	umax = xPixels = x;
	vmax = yPixels = y;
	super.setSize(x, y);
    }

    protected double xmin_ = -1;
    protected double xmax_ = 1;
    protected double ymin_ = -1;
    protected double ymax_ = 1;

    /**
     * Sets up a mapping from physical units to pixel units.
     *
     * @param xmin the left edge of the physical window
     * @param xmax the right edge of the physical window
     * @param ymin the bottom edge of the physical window
     * @param ymax the top edge of the physical window
     */
    public void setWindow (double xmin, double xmax,
			   double ymin, double ymax) {
	this.xmin_ = xmin;
	this.xmax_ = xmax;
	this.ymin_ = ymin;
	this.ymax_ = ymax;
    }

    protected int umin = 0;
    protected int umax = xPixels;
    protected int vmin = 0;
    protected int vmax = yPixels;

    /**
     * Computes the x pixel value of a point in the world window.
     *
     * @param x the x coordinate of the point in the world window
     * @param y the y coordinate of the point in the world window
     * @return the x pixel value of this point on the plot
     */
    public int xPixel (double x, double y) {
	double u = umin + (x - xmin_) / (xmax_ - xmin_) * (umax - umin);
	return (int) Math.round(u);
    }

    public double xWorld (int x) {
	double wx = x - umin;
	wx /= umax - umin;
	return xmin_ + wx * (xmax_ - xmin_);
    }

    public double yWorld (int y) {
	double wy = y - vmin;
	wy /= vmax - vmin;
	return ymax_ - wy * (ymax_ - ymin_);
    }

    /**
     * Computes the y pixel value of a point in the world window.
     *
     * @param x the x coordinate of the point in the world window
     * @param y the y coordinate of the point in the world window
     * @return the y pixel value of this point on the plot
     */
    public int yPixel (double x, double y) {
	double v = vmin + (ymax_ - y) / (ymax_ - ymin_) * (vmax - vmin);
	return (int) Math.round(v);
    }

    /**
     * Computes the (x, y) pixel values of a point in the world window.
     *
     * @param x the x coordinate of the point in the world window
     * @param y the y coordinate of the point in the world window
     * @return the (x, y) pixel value of this point on the plot
     */
    public Point getPoint (double x, double y) {
	return new Point(xPixel(x, y), yPixel(x, y));
    }

    /**
     * Computes the separation in pixels of two points in the world window.
     *
     * @param dx the x coordinate difference of the two points in the world
     * @param dy the y coordinate difference of the two points in the world
     * @return the separation in pixels
     */
    public Dimension getDimension (double dx, double dy) {
	return new Dimension(xPixel(dx, dy) - xPixel(0, 0),
			     yPixel(0, 0) - yPixel(dx, dy));
    }

    protected Image offScreen, offScreen2;
    protected Graphics osg, osg2;
    protected boolean doOverlay;
    
    public void plotOverlay () {
	doOverlay = true;
	if (offScreen2 == null) {
	    offScreen2 = createImage(getSize().width,
				    getSize().height);
	    osg2 = offScreen2.getGraphics();
	}
	osg2.setPaintMode();
	osg2.drawImage(offScreen, 0, 0, null);
	osg2.setXORMode(getBackground());
	Graphics swap = osg;
	osg = osg2;
	osg2 = swap;
    }

    /**
     * Classes which inherit from Plot must define this method.
     * This method is called when the plot is repainted.
     * All instructions which draw in the plot should be placed in this method.
     */
    public abstract void paint ();

    /**
     * Clears the plot.
     */
    public void clear () {
	osg.clearRect(0, 0, getSize().width, getSize().height);
    }

    public void plotPoint (double x, double y) {
	int ix = xPixel(x, y);
	int iy = yPixel(x, y);
	osg.drawLine(ix, iy, ix, iy);
    }

    public void plotLine (double x1, double y1, double x2, double y2) {
	int ix1 = xPixel(x1, y1);
	int iy1 = yPixel(x1, y1);
	int ix2 = xPixel(x2, y2);
	int iy2 = yPixel(x2, y2);
	osg.drawLine(ix1, iy1, ix2, iy2);
    }

    public void plotString (String s, double x, double y) {
	int ix = xPixel(x, y);
	int iy = yPixel(x, y);
	osg.drawString(s, ix, iy);
    }

    public void plotStringCenter (String s, double x, double y) {
	int ix = xPixel(x, y);
	int iy = yPixel(x, y);
	FontMetrics fm = osg.getFontMetrics();
	int sw = fm.stringWidth(s);
	ix -= sw / 2;
	osg.drawString(s, ix, iy);
    }

    public void plotStringLeft (String s, double x, double y) {
	int ix = xPixel(x, y);
	int iy = yPixel(x, y);
	FontMetrics fm = osg.getFontMetrics();
	int sw = fm.stringWidth(s);
	ix -= sw;
	osg.drawString(s, ix, iy);
    }

    public void boxArea (double x1, double x2, double y1, double y2) {
	double x = Math.min(x1, x2);
	double y = Math.max(y1, y2);
	double dx = Math.abs(x1 - x2);
	double dy = Math.abs(y1 - y2);
	int ix = xPixel(x, y);
	int iy = yPixel(x, y);
	Dimension d = getDimension(dx, dy);
	osg.fillRect(ix, iy, d.width, d.height);
    }

    public void plotBox (double x1, double x2, double y1, double y2) {
	double x = Math.min(x1, x2);
	double y = Math.max(y1, y2);
	double dx = Math.abs(x1 - x2);
	double dy = Math.abs(y1 - y2);
	int ix = xPixel(x, y);
	int iy = yPixel(x, y);
	Dimension d = getDimension(dx, dy);
	osg.drawRect(ix, iy, d.width, d.height);
    }

    public void floodCircle (double x1, double x2, double y1, double y2) {
	double x = Math.min(x1, x2);
	double y = Math.max(y1, y2);
	double dx = Math.abs(x1 - x2);
	double dy = Math.abs(y1 - y2);
	int ix = xPixel(x, y);
	int iy = yPixel(x, y);
	Dimension d = getDimension(dx, dy);
	osg.fillOval(ix, iy, d.width, d.height);
    }
    
    public void circle (double x1, double x2, double y1, double y2) {
	double x = Math.min(x1, x2);
	double y = Math.max(y1, y2);
	double dx = Math.abs(x1 - x2);
	double dy = Math.abs(y1 - y2);
	int ix = xPixel(x, y);
	int iy = yPixel(x, y);
	Dimension d = getDimension(dx, dy);
	osg.drawOval(ix, iy, d.width, d.height);
    }

    public void boxClear (double x1, double x2, double y1, double y2) {
	double x = Math.min(x1, x2);
	double y = Math.max(y1, y2);
	double dx = Math.abs(x1 - x2);
	double dy = Math.abs(y1 - y2);
	int ix = xPixel(x, y);
	int iy = yPixel(x, y);
	Dimension d = getDimension(dx, dy);
	osg.clearRect(ix, iy, d.width, d.height);
    }

    public void update (Graphics g) {
	paint(g);
    }

    public void paint (Graphics g) {
	if (offScreen == null) {
	    offScreen = createImage(getSize().width,
				    getSize().height);
	    osg = offScreen.getGraphics();
	}
	paint();
	if (!doOverlay) {
	    g.drawImage(offScreen, 0, 0, null);
	} else {
	    g.drawImage(offScreen2, 0, 0, null);
	    Graphics swap = osg2;
	    osg2 = osg;
	    osg = swap;
	    doOverlay = false;
	}
    }

    public void drawAxes (double xmin, double xmax,
			  double ymin, double ymax) {
	int ntick = 10;
	double dx = (xmax - xmin) / ntick;
	double dy = (ymax - ymin) / ntick;
//	setWindow(xmin - dx, xmax + dx, ymin - dy, ymax + dy);
	this.xmin_ = xmin - dx;
	this.xmax_ = xmax + dx;
	this.ymin_ = ymin - dy;
	this.ymax_ = ymax + dy;

	double x0 = Math.max(0, xmin);
	double y0 = Math.max(0, ymin);
	if (ymin * ymax < 0)
	    y0 = 0;
	else
	    y0 = ymin;
	plotLine(xmin, y0, xmax, y0);
	plotLine(x0, ymin, x0, ymax);
	double tx = 0.1 * dy;
	double ty = 0.1 * dx;
	for (int itick = 0; itick <= ntick; itick++) {
	    double x = xmin + itick * dx;
	    plotLine(x, y0 - tx, x, y0 + tx);
	    double y = ymin + itick * dy;
	    plotLine(x0 - ty, y, x0 + ty, y);
	}
    }

    public void setColor (Color c) {
	osg.setColor(c);
    }

    Hashtable colorTable;
    String[] awtColorString
	= {"black", "blue", "cyan", "darkGray", "gray", "green", "lightGray",
	   "magenta", "orange", "pink", "red", "white", "yellow"};
    Color[] awtColor
	= {Color.black, Color.blue, Color.cyan, Color.darkGray, Color.gray,
	   Color.green, Color.lightGray, Color.magenta, Color.orange,
	   Color.pink, Color.red, Color.white, Color.yellow};

    private void populateColorTable () {
	for (int i = 0; i < awtColor.length; i++)
	    colorTable.put(awtColorString[i], awtColor[i]);
    }

    public void setColor (String s) {
	if (colorTable == null) {
	    colorTable = new Hashtable();
	    populateColorTable();
	}
	Color c = (Color) colorTable.get(s);
	if (c != null)
	    osg.setColor(c);
	else
	    System.err.println("comphys.graphics.Plot: no such color " + s);
    }

    public void setBackground (String s) {
	if (colorTable == null) {
	    colorTable = new Hashtable();
	    populateColorTable();
	}
	Color c = (Color) colorTable.get(s);
	if (c != null)
	    setBackground(c);
	else
	    System.err.println("comphys.graphics.Plot: no such color " + s);
    }

    public void setXORMode () {
	osg.setXORMode(getBackground());
    }

    public void setPaintMode () {
	osg.setPaintMode();
    }

}
