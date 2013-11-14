package g3dsys.images;


/**
 *	Superclass for drawable shapes.
 *	The subclasses store all location information as molecule space coordinates
 *  (in Angstroms). These are converted to pixel coordinates based on
 *  current rotation, translation, zoom, and perspective settings by the
 *  CoordMapper whenever the Figures need to draw themselves.
 *  
 *  The drawing methods are intended to be called by FigureManager.
 */

public abstract class Figure {
	
	protected int _id; //serial id for figure tracking and removal
	
    boolean drawme = true;
    
	protected g3dsys.control.G3DSys _gsys;
	
	Figure(g3dsys.control.G3DSys g) {
		_gsys = g;
        _id = -1;
	}
	
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
