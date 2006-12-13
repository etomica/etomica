package g3dsys.images;

import javax.vecmath.Point3f;

/**
 *	Superclass for drawable shapes.
 *	All location information is stored as molecule space coordinates
 *  (in Angstroms). These are converted to pixel coordinates based on
 *  current rotation, translation, zoom, and perspective settings by the
 *  CoordMapper whenever the Figures need to draw themselves.
 *  
 *  For normal use the Figure constructors should not be called directly,
 *  as G3DSys will do this itself. The drawing methods are intended to be
 *  called by FigureManager.
 */

public abstract class Figure {
	
	protected short _c; //color
	protected Point3f _p; //location in molecule space
	protected float _d; //diameter of the figure, when applicable
	protected int _id; //serial id for figure tracking and removal
	
    boolean drawme = true;
    
	protected g3dsys.control.G3DSys _gsys;
	
	Figure(g3dsys.control.G3DSys g, short c) {
		_gsys = g;
		_c = c;
		_p = new Point3f(0,0,0); //default location is origin
	}
	
	/** @return the color of the figure */
	public short getColor() { return _c; }
	/** set the color of the figure */
	public void setColor(short c) { _c = c; }
	/** @return the molspace position of the figure */
	public Point3f getPoint() { return _p; }
	/** set the molspace location of the figure */
	public void setPoint(Point3f p) { _p = p; }
	/** @return the x molspace position of the figure */
	public float getX() { return _p.x; }
	/** set the x molspace position of the figure */
	public void setX(float x) { _p.x = x; }
	/** @return the y molspace position of the figure */
	public float getY() { return _p.y; }
	/** set the y molspace position of the figure */
	public void setY(float y) { _p.y = y; }
	/** @return the z molspace position of the figure */
	public float getZ() { return _p.z; }
	/** set the z molspace position of the figure */
	public void setZ(float z) { _p.z = z; }
	/** @return the 'size' (largest dimension) of the figure in Angstroms */
	public abstract float getD();
	/** set the molspace diameter of the figure, when applicable */
	public void setD(float d) { _d = d; }
	
	/**Get the ID of the figure
	 * @return the ID of the figure
	 */
	public int getID() { return _id; }
	/**
	 * Set the ID of the figure. Should only be used by the FigureManager
	 * @param id the ID to set
	 */
	public void setID(int id) { _id = id; }
	
	/**
	 * All Figures know how to draw themselves via override.
	 */
	abstract public void draw();
	
    public void setDrawable(boolean drawable) { drawme = drawable; }
    
}
