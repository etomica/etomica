/*
 * Created on Jun 15, 2005
 */
package etomica.models.hexane;

import etomica.space.Space;
import etomica.space.BoundaryDeformablePeriodic;
import etomica.space3d.Vector3D;

/**
 * Class that defines the boundary for a cell based on Dr. Monson's data.
 * 
 * @author nancycribbin
 *
 */
public class BoundaryHexane extends BoundaryDeformablePeriodic {
    
    public BoundaryHexane(Space space, boolean[] periodicity){
        super(space, periodicity, 
                new Vector3D[] {new Vector3D(5.79141400642806, 6.37543175848628, 12.16404681464840),
                new Vector3D(0.23085087646072, 0.25413033978032, 0.48486964761179),
                new Vector3D(0.13790821930218, 1.21905982551516, 2.20715141449381)});
        makePeriodicity(space.D());
    }
    
    
    
    /**
     * Periodic in all directions; same as the superclass's private method.
     * Same as th
     * @param D
     * @return
     */
	private static boolean[] makePeriodicity(int D) {
		boolean[] isPeriodic = new boolean[D];
		for (int i = 0; i < D; i++) {
			isPeriodic[i] = true;
		}
		return isPeriodic;
	}
}
