package g3dsys.images;

import g3dsys.control.G3DSys;

import javax.vecmath.Point3f;
import javax.vecmath.Point3i;

import org.jmol.g3d.Graphics3D;


/**
 * Class for drawing a bounding box around the figures in molecule space
 */

public class Box extends Figure {
	
	// eight points to define the box
	private Point3f LUT,LDT,RUT,RDT; //leftUpTop, leftDownTop, etc.
	private Point3f LUB,LDB,RUB,RDB; //leftUpBottom, leftDownBottom, etc.
    private final Point3i t1, t2;

	public Box(G3DSys g, short c) {
		super(g, c);
		
		resize();
		
		_p = new Point3f(
				(_gsys.getMinX() + _gsys.getMaxX())/2,
				(_gsys.getMinY() + _gsys.getMaxY())/2,
				(_gsys.getMinZ() + _gsys.getMaxZ())/2);
        
        t1 = new Point3i();
        t2 = new Point3i();
		
		g.setCenterOfRotation(_p); //rotate about center of box
	}

	private void resize() {
		float minX = _gsys.getMinX(); float maxX = _gsys.getMaxX();
		float minY = _gsys.getMinY(); float maxY = _gsys.getMaxY();
		float minZ = _gsys.getMinZ(); float maxZ = _gsys.getMaxZ();

		LUT = new Point3f(minX,minY,minZ); LUB = new Point3f(minX,minY,maxZ);
		LDT = new Point3f(minX,maxY,minZ); LDB = new Point3f(minX,maxY,maxZ);
		RUT = new Point3f(maxX,minY,minZ); RUB = new Point3f(maxX,minY,maxZ);
		RDT = new Point3f(maxX,maxY,minZ); RDB = new Point3f(maxX,maxY,maxZ);
	}

	public void draw() {
		resize(); //must resize in case the model has changed
		Graphics3D g3d = _gsys.getG3D();
		
		//four top to bottom lines
        _gsys.screenSpace(LUT, t1); _gsys.screenSpace(LUB, t2);
		g3d.drawLine(_c, t1, t2);
        _gsys.screenSpace(LDT, t1); _gsys.screenSpace(LDB, t2);
        g3d.drawLine(_c, t1, t2);
        _gsys.screenSpace(RUT, t1); _gsys.screenSpace(RUB, t2);
        g3d.drawLine(_c, t1, t2);
        _gsys.screenSpace(RDT, t1); _gsys.screenSpace(RDB, t2);
        g3d.drawLine(_c, t1, t2);
		
		//four left to right lines
        _gsys.screenSpace(LUT, t1); _gsys.screenSpace(RUT, t2);
        g3d.drawLine(_c, t1, t2);
        _gsys.screenSpace(LDT, t1); _gsys.screenSpace(RDT, t2);
        g3d.drawLine(_c, t1, t2);
        _gsys.screenSpace(LUB, t1); _gsys.screenSpace(RUB, t2);
        g3d.drawLine(_c, t1, t2);
        _gsys.screenSpace(LDB, t1); _gsys.screenSpace(RDB, t2);
        g3d.drawLine(_c, t1, t2);
		
		//four up to down lines
        _gsys.screenSpace(LUT, t1); _gsys.screenSpace(LDT, t2);
        g3d.drawLine(_c, t1, t2);
        _gsys.screenSpace(RUT, t1); _gsys.screenSpace(RDT, t2);
        g3d.drawLine(_c, t1, t2);
        _gsys.screenSpace(LUB, t1); _gsys.screenSpace(LDB, t2);
        g3d.drawLine(_c, t1, t2);
        _gsys.screenSpace(RUB, t1); _gsys.screenSpace(RDB, t2);
        g3d.drawLine(_c, t1, t2);
	}
	
	public float getD() {
		return 0;
	}

}
