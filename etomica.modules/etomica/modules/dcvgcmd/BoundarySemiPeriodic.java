package etomica.modules.dcvgcmd;

import etomica.Space;
import etomica.space.BoundaryPeriodicSquare;
import etomica.space.Coordinate;
import etomica.space3d.Vector3D;

/**
 * @author kofke
 *
 * To change this generated comment edit the template variable "typecomment":
 * Window>Preferences>Java>Templates.
 * To enable and disable the creation of type comments go to
 * Window>Preferences>Java>Code Generation.
 */
public class BoundarySemiPeriodic extends BoundaryPeriodicSquare {
	
	public BoundarySemiPeriodic(Space space){
		super(space);
	}
	
	public void nearestImage(Vector3D dr) {
  //      dr.x -= dimensions.x*((dr.x > 0.0) ? Math.floor(dr.x/dimensions.x + 0.5) : Math.ceil(dr.x/dimensions.x - 0.5));
  //      dr.y -= dimensions.y*((dr.y > 0.0) ? Math.floor(dr.y/dimensions.y + 0.5) : Math.ceil(dr.y/dimensions.y - 0.5));
  //      dr.z -= dimensions.z*((dr.z > 0.0) ? Math.floor(dr.z/dimensions.z + 0.5) : Math.ceil(dr.z/dimensions.z - 0.5));
  //      final double dimxHalf = 0.5*dimensions.x;
  //      final double dimyHalf = 0.5*dimensions.y;
  //      final double dimzHalf = 0.5*dimensions.z;
	/*    while(dr.x > +dimensionsHalf.x) dr.x -= dimensions.x;
		while(dr.x < -dimensionsHalf.x) dr.x += dimensions.x;*/
		while(dr.x(1) > +dimensionsHalf.x(1)) dr.setX(1,dr.x(1)-dimensions.x(1));
		while(dr.x(1) < -dimensionsHalf.x(1)) dr.setX(1,dr.x(1)+dimensions.x(1));
		while(dr.x(0) > +dimensionsHalf.x(0)) dr.setX(0,dr.x(0)-dimensions.x(0));
		while(dr.x(0) < -dimensionsHalf.x(0)) dr.setX(0,dr.x(0)+dimensions.x(0));
		//System.out.println("dimesions = "+dimensions);
	//    System.out.print(dr.x+"  ");dr.x %= dimensionsHalf.x; System.out.println(dr.x);
	//    dr.x = ((dr.x + dimensions.x) % dimensions.x) - dimensionsHalf.x;
	//    dr.y = ((dr.y + dimensions.y) % dimensions.y) - dimensionsHalf.y;
	//    dr.z = ((dr.z + dimensions.z) % dimensions.z) - dimensionsHalf.z;
	}

	public boolean centralImage(Coordinate c) {
		Vector3D r = (Vector3D) c.position();
		if(r.x(0)<0.0) r.setX(0,r.x(0)+dimensions.x(0));
		if(r.x(0)>dimensions.x(0)) r.setX(0,r.x(0)-dimensions.x(0));
		if(r.x(1)<0.0) r.setX(1,r.x(1)+dimensions.x(1));
		if(r.x(1)>dimensions.x(1)) r.setX(1,r.x(1)-dimensions.x(1));
		return true;
	}

}
