package etomica.modules.clustergenerator;

import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.util.Hashtable;

/**
 * 
 * Used from comphys.graphic.plot
 * A class which makes plotting a little more convenient than using
 * Java's Graphics methods.
 * 
 */
public abstract class Plot extends Component {

    /**
     * Constructor sets default size of the plot
     */
    public Plot () {
        this(Color.white,400,400);
    }
    
    public Plot(Color background, int xSize, int ySize) {
        setBackground(background);
        setSize(xSize,ySize);
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


    /**
     * Sets up a mapping from physical units to pixel units.
     *
     * @param xmin the left edge of the physical window
     * @param xmax the right edge of the physical window
     * @param ymin the bottom edge of the physical window
     * @param ymax the top edge of the physical window
     */
    public void setWindow(double xmin, double xmax, double ymin, double ymax) {
        this.xmin_ = xmin;
        this.xmax_ = xmax;
        this.ymin_ = ymin;
        this.ymax_ = ymax;
    }

    /**
     * Computes the x pixel value of a point in the world window.
     *
     * @param x the x coordinate of the point in the world window
     * @param y the y coordinate of the point in the world window
     * @return the x pixel value of this point on the plot
     */
    public int xPixel(double x) {
        double u = umin + (x - xmin_) / (xmax_ - xmin_) * (umax - umin);
        return (int) Math.round(u);
    }

    /**
     * Computes the y pixel value of a point in the world window.
     *
     * @param x the x coordinate of the point in the world window
     * @param y the y coordinate of the point in the world window
     * @return the y pixel value of this point on the plot
     */
    public int yPixel(double y) {
        double v = vmin + (ymax_ - y) / (ymax_ - ymin_) * (vmax - vmin);
        return (int) Math.round(v);
    }

    /**
     * Computes the separation in pixels of two points in the world window.
     *
     * @param dx the x coordinate difference of the two points in the world
     * @param dy the y coordinate difference of the two points in the world
     * @return the separation in pixels
     */
    public Dimension getDimension(double dx, double dy) {
        return new Dimension(xPixel(dx) - xPixel(0), yPixel(0) - yPixel(dy));
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
    public void clear() {
        osg.clearRect(0, 0, getSize().width, getSize().height);
    }

    public void plotLine(double x1, double y1, double x2, double y2) {
        int ix1 = xPixel(x1);
        int iy1 = yPixel(y1);
        int ix2 = xPixel(x2);
        int iy2 = yPixel(y2);
        osg.drawLine(ix1, iy1, ix2, iy2);
    }

    public void plotString(String s, double x, double y) {
        int ix = xPixel(x);
        int iy = yPixel(y);
        osg.drawString(s, ix, iy);
    }

    public void plotStringCenter(String s, double x, double y) {
        int ix = xPixel(x);
        int iy = yPixel(y);
        FontMetrics fm = osg.getFontMetrics();
        int sw = fm.stringWidth(s);
        ix -= sw / 2;
        osg.drawString(s, ix, iy);
    }

    public void plotBox(double x1, double x2, double y1, double y2) {
        double x = Math.min(x1, x2);
        double y = Math.max(y1, y2);
        double dx = Math.abs(x1 - x2);
        double dy = Math.abs(y1 - y2);
        int ix = xPixel(x);
        int iy = yPixel(y);
        Dimension d = getDimension(dx, dy);
        osg.drawRect(ix, iy, d.width, d.height);
    }

    public void floodInside (double x1, double x2, double y1, double y2) {
    	double x = Math.min(x1, x2);
    	double y = Math.max(y1, y2);
    	double dx = Math.abs(x1 - x2)-0.1;
    	double dy = Math.abs(y1 - y2)-0.1;
    	int ix = xPixel(x);
    	int iy = yPixel(y);
    	Dimension d = getDimension(dx, dy);
    	osg.fillOval(ix, iy, d.width, d.height);
        }

    public void floodCircle(double x1, double x2, double y1, double y2) {
        double x = Math.min(x1, x2);
        double y = Math.max(y1, y2);
        double dx = Math.abs(x1 - x2);
        double dy = Math.abs(y1 - y2);
        int ix = xPixel(x);
        int iy = yPixel(y);
        Dimension d = getDimension(dx, dy);
        osg.fillOval(ix, iy, d.width, d.height);
    }
    
    public void circle(double x1, double x2, double y1, double y2) {
        double x = Math.min(x1, x2);
        double y = Math.max(y1, y2);
        double dx = Math.abs(x1 - x2);
        double dy = Math.abs(y1 - y2);
        int ix = xPixel(x);
        int iy = yPixel(y);
        Dimension d = getDimension(dx, dy);
        osg.drawOval(ix, iy, d.width, d.height);
    }

    public void update(Graphics g) {
        paint(g);
    }

    public void paint (Graphics g) {
        osg = g;
        paint();
    }

    private void populateColorTable() {
        for (int i = 0; i < awtColor.length; i++)
            colorTable.put(awtColorString[i], awtColor[i]);
    }

    public void setColor(String s) {
        if (colorTable == null) {
            colorTable = new Hashtable();
            populateColorTable();
        }
        Color c = (Color) colorTable.get(s);
        if (c != null)
            osg.setColor(c);
        else
            System.err.println("Plot: no such color " + s);
    }

    protected int xPixels;
    protected int yPixels;
    
    protected double xmin_ = -1;
    protected double xmax_ = 1;
    protected double ymin_ = -1;
    protected double ymax_ = 1;

    protected int umin;
    protected int umax;
    protected int vmin;
    protected int vmax;
    
    protected Graphics osg;

    Hashtable colorTable;
    String[] awtColorString = { "black", "blue", "cyan", "darkGray", "gray",
            "green", "lightGray", "magenta", "orange", "pink", "red", "white", "yellow"};
    Color[] awtColor = { Color.black, Color.blue, Color.cyan, Color.darkGray,
            Color.gray, Color.green, Color.lightGray, Color.magenta,
            Color.orange, Color.pink, Color.red, Color.white, Color.yellow};

}
