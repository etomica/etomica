package etomica.math.geometry.coordinate;

import etomica.space.IVector;
import etomica.space1d.Vector1D;
import etomica.space2d.Vector2D;
import etomica.space3d.Vector3D;

public class CoordinateConverter {

    public static void toSpherical(IVector v, double[] result) {
        switch (v.getD()) {
            case 1:
                toSpherical((Vector1D)v, result);
                break;
            case 2:
                toSpherical((Vector2D)v, result);
                break;
            case 3:
                toSpherical((Vector3D)v, result);
                break;
            default:
                throw new IllegalArgumentException("I'm impressed by your ability to make a "+v.getD()+" dimension vector.  Please try 1-3 next time");
        }
    }
    
    public static void toSpherical(Vector1D v, double[] result) {
        result[0] = v.x(0);
    }
    
    public static void toSpherical(Vector2D v, double[] result) {
        result[0] = Math.sqrt(v.squared());
        result[1] = Math.acos(v.x(1) / result[0]); //theta
    }
    
    public static void toSpherical(Vector3D v, double[] result) {
        result[0] = Math.sqrt(v.squared());
        result[1] = Math.acos(v.x(2) / result[0]); //theta
        result[2] = Math.atan2(v.x(2), v.x(1)); //phi
    }

}
