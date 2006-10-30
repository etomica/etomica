package g3dsys.images;

import g3dsys.control.G3DSys;

import javax.vecmath.Point3f;

import org.jmol.g3d.Graphics3D;


/**
 * Class for drawing a bounding box around the figures in molecule space
 */

public class Box extends Figure {
	
	// eight points to define the box
	private Point3f LUT,LDT,RUT,RDT; //leftUpTop, leftDownTop, etc.
	private Point3f LUB,LDB,RUB,RDB; //leftUpBottom, leftDownBottom, etc.

	public Box(G3DSys g, short c) {
		super(g, c);
		
		resize();
		
		_p = new Point3f(
				(_gsys.getMinX() + _gsys.getMaxX())/2,
				(_gsys.getMinY() + _gsys.getMaxY())/2,
				(_gsys.getMinZ() + _gsys.getMaxZ())/2);
		
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
		
		g3d.fillSphereCentered(org.jmol.g3d.Graphics3D.GOLD,
				10, _gsys.screenSpace(_p));
		
		//four top to bottom lines
		g3d.drawLine(_c, _gsys.screenSpace(LUT), _gsys.screenSpace(LUB));
		g3d.drawLine(_c, _gsys.screenSpace(LDT), _gsys.screenSpace(LDB));
		g3d.drawLine(_c, _gsys.screenSpace(RUT), _gsys.screenSpace(RUB));
		g3d.drawLine(_c, _gsys.screenSpace(RDT), _gsys.screenSpace(RDB));
		
		//four left to right lines
		g3d.drawLine(_c, _gsys.screenSpace(LUT), _gsys.screenSpace(RUT));
		g3d.drawLine(_c, _gsys.screenSpace(LDT), _gsys.screenSpace(RDT));
		g3d.drawLine(_c, _gsys.screenSpace(LUB), _gsys.screenSpace(RUB));
		g3d.drawLine(_c, _gsys.screenSpace(LDB), _gsys.screenSpace(RDB));
		
		//four up to down lines
		g3d.drawLine(_c, _gsys.screenSpace(LUT), _gsys.screenSpace(LDT));
		g3d.drawLine(_c, _gsys.screenSpace(RUT), _gsys.screenSpace(RDT));
		g3d.drawLine(_c, _gsys.screenSpace(LUB), _gsys.screenSpace(LDB));
		g3d.drawLine(_c, _gsys.screenSpace(RUB), _gsys.screenSpace(RDB));		
	}
	
	public float getD() {
		return 0;
	}

}
