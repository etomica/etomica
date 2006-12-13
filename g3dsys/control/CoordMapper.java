package g3dsys.control;

import javax.vecmath.Matrix3f;
import javax.vecmath.Point3f;
import javax.vecmath.Point3i;

/**
 *	Class that maps from molecule space to the pixel space
 *	of g3d drawing functions. Handles scaling for perspective, resize,
 *  and zoom as well.
 */

class CoordMapper {
	
	// pixels per angstrom, determined by molecule size, zoom level, etc.
	private float ppa;

	// zoom level
	private float zoomLevel;
	
	// rotation matrix
	private Matrix3f rotation;
	private Matrix3f tmp = new Matrix3f();
	
	// comm channel
	private G3DSys master;
	
	// translation offsets
	private int xlateHoriz, xlateVert;
	
	// center of rotation in Angstroms; origin by default
	private Point3f cor = new Point3f(0,0,0);
    private Point3f tempP = new Point3f();
	
	//perspective flag and scaling factor
	private boolean PERSPECTIVE = true;
	private float PERSPECTIVE_SCALE = .03f;
	
	public CoordMapper(G3DSys m) {
		ppa = 1f;
		rotation = new Matrix3f(1,0,0,0,1,0,0,0,1);
		xlateHoriz = xlateVert = 0;
		zoomLevel = 100f;
		master = m;
	}

	/* *****************************************************************
	 * Rotation methods
	 * *****************************************************************/

	/** Rotate to the home position; remove all rotations and translations */
	public void rotateToHome() {
		rotation.setIdentity();
		xlateHoriz = xlateVert = 0;
	}
	
	/**Set a specific x axis rotation
	 * @param degrees rotation amount to set, in degrees
	 */
	public void rotateToX(float degrees) { rotation.rotX(degToRad(degrees)); }
	
	/**Set a specific y axis rotation
	 * @param degrees rotation amount to set, in degrees
	 */
	public void rotateToY(float degrees) { rotation.rotY(degToRad(degrees)); }
	
	/**Set a specific z axis rotation
	 * @param degrees rotation amount to set, in degrees
	 */
	public void rotateToZ(float degrees) { rotation.rotY(degToRad(degrees)); }
	
	/**Rotate about the x axis relative to current value
	 * @param degrees rotation amount to add, in degrees
	 */
	public void rotateByX(float degrees) {		
		tmp.rotX( degToRad(degrees) );
		rotation.mul(tmp,rotation);
	}
	
	/**Rotate about the y axis relative to current value
	 * @param degrees rotation amount to add, in degrees
	 */
	public void rotateByY(float degrees) {
		tmp.rotY( degToRad(degrees) );
		rotation.mul(tmp,rotation);
	}
	
	/**Rotate about the z axis relative to current value
	 * @param degrees rotation amount to add, in degrees
	 */
	public void rotateByZ(float degrees) {
		tmp.rotZ( degToRad(degrees) );
		rotation.mul(tmp,rotation);
	}
	/**Convert degrees to radians
	 * @param degrees the degree value to convert
	 * @return the radian value
	 */
	private float degToRad(float degrees) {
		return degrees * (float)java.lang.Math.PI/180;
	}

	/**Set the point around which rotation will take place
	 * @param p the point around which to rotate
	 */
	public void setCenterOfRotation(Point3f p) { cor = p; }
	
	/**Get the point around which rotation takes place
	 * @return the center of rotation
	 */
	public Point3f getCenterOfRotation() { return cor; }

	
	
	/* *****************************************************************
	 * Display and scaling methods
	 * *****************************************************************/
	
	/**Convert from Angstroms to pixels
	 * @param angstroms number to convert
	 * @return the number of pixels needed
	 */
	public int angToPixel(float angstroms) {
		return (int)(angstroms * ppa * (zoomLevel/100f));
	}
	
	/**Gives G3D screen coordinate to paint to for a given molecule coordinate
	 * @param molSpace the molecule space coordinate in Angstroms
	 * @return the G3D pixel space coordinate after conversion and rotation
	 */
	public void screenSpace(Point3f molSpace, Point3i screenSpace) {
		//translate and offset for the center of rotation
        tempP.x = molSpace.x - cor.x;
        tempP.y = molSpace.y - cor.y;
        tempP.z = molSpace.z - cor.z;
		
		//rotate and remove rotation offsets
        rotation.transform(tempP);
		tempP.x += cor.x;
		tempP.y += cor.y;
		tempP.z += cor.z;
		
		//need better centering offsets
		//-- these seem to work well now
		screenSpace.x = xlateHoriz +
			angToPixel(tempP.x + java.lang.Math.abs(master.getMinX()));
		screenSpace.y = xlateVert + 
			angToPixel(tempP.y + java.lang.Math.abs(master.getMinY()));
		screenSpace.z = angToPixel(tempP.z) + (int) Short.MAX_VALUE/2;

	}

	/**Performs perspective calculation, if enabled
	 * @param z pixel depth of the object
	 * @param diameter internal diameter in Angstroms
	 * @return drawing diameter in pixels
	 */
	public int perspective(int z, float diameter) {
		int frontSide = Short.MAX_VALUE/2 - angToPixel(master.getAngstromDepth()/2);
		//System.out.println("slab starts at depth "+frontSide+", we are "+z);
		if(PERSPECTIVE)
			return (angToPixel(diameter)-(int)((z - frontSide)*PERSPECTIVE_SCALE));
		return angToPixel(diameter);
	}

	/**Set the pixel to Angstrom ratio
	 * @param f the ppa ratio to use
	 */
	public void setPPA(float f) { ppa = f; }
	
	/**Get the pixel to Angstrom ratio
	 * @return the ppa ratio
	 */
	public float getPPA() { return ppa; }
	
	/**
	 * Recalculates the PPA based on model size in Angstroms and the available
	 * window size in pixels.
	 */
	public void recalcPPA() {
		//TODO: zoom level should impact this, instead of diameter as it does now

		//System.out.print("window size "+master.getPixelWidth()+"x"+master.getPixelHeight());
		//System.out.println("; mol size "+master.getAngstromWidth()+"x"+master.getAngstromHeight());
		int smallDim =
			java.lang.Math.min(master.getPixelWidth(), master.getPixelHeight());
		float bigDim = java.lang.Math.max(
				master.getAngstromWidth(), master.getAngstromHeight());
		ppa = smallDim/bigDim;
	}

	
	
	
	/* *****************************************************************
	 * Zoom
	 * *****************************************************************/
	/** Increase the zoom level i% */
	public void zoomUp(int i) { zoomLevel += i; } //yes, we can 'increase' by -x%
	/** Decrease the zoom level i% */
	public void zoomDown(int i) { zoomLevel -= i; }//yes, we can 'decrease' by -x%
	
	
	
	
	/* *****************************************************************
	 * Translation
	 * *****************************************************************/
	/** Translate horizontally x pixels */
	public void xlateX(int x) { xlateHoriz += x; }
	/** Translate vertically y pixels */
	public void xlateY(int y) { xlateVert += y; }
	
}
