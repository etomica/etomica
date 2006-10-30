package g3dsys.images;

import g3dsys.control.G3DSys;

import javax.vecmath.Point3f;
import javax.vecmath.Point3i;


import org.jmol.g3d.Graphics3D;

/**
 * Class for drawing cylinders (bonds)
 */

public class Cylinder extends Figure {
	
	private Point3f p1; //center of first endcap
	private Point3f p2; //center of second endcap

	public Cylinder(G3DSys g, short c, float x1, float y1, float z1,
			float x2, float y2, float z2) {
		this(g, c, new Point3f(x1,y1,z1), new Point3f(x2,y2,z2));
	}
	public Cylinder(G3DSys g, short c, Point3f p1, Point3f p2) {
		super(g,c);
		this.p1 = p1;
		this.p2 = p2;
		_p = new Point3f( (p1.x+p2.x)/2, (p1.y+p2.y)/2, (p1.z+p2.z)/2 );
	}

	public void draw() {
		Point3i p1i = _gsys.screenSpace(p1);
		Point3i p2i = _gsys.screenSpace(p2);
		_gsys.getG3D().fillCylinder(_c, Graphics3D.ENDCAPS_FLAT, dist(p1i,p2i)/2, p1i, p2i);
	}
	
	/**
	 * @return distance between points defining cylinder, within one pixel
	 */
	private int dist(Point3i p1i, Point3i p2i) {
		return (int) java.lang.Math.sqrt((p2i.x-p1i.x)*(p2i.x-p1i.x) + (p2i.y-p1i.y)*(p2i.y-p1i.y) + (p2i.z-p1i.z)*(p2i.z-p1i.z)); 
	}
	/**
	 * @return distance between points defining cylinder in molecule space
	 */
	private float dist(Point3f p1f, Point3f p2f) {
		return (float) java.lang.Math.sqrt((p2f.x-p1f.x)*(p2f.x-p1f.x) + (p2f.y-p1f.y)*(p2f.y-p1f.y) + (p2f.z-p1f.z)*(p2f.z-p1f.z));
	}
	
	/**Get the 'diameter' of the cylinder in Angstroms
	 * @return the 'diameter' of the cylinder; the distance between endpoints
	 */
	public float getD() {
		return (float)dist(p1,p2);
	}

}
