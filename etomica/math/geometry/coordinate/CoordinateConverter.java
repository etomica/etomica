package etomica.math.geometry.coordinate;

import etomica.api.IVector;

public class CoordinateConverter {

    public static void toSpherical(IVector v, double[] result) {
        switch (v.getD()) {
            case 1:
                result[0] = v.x(0);
                break;
            case 2:
                result[0] = Math.sqrt(v.squared());
                result[1] = Math.acos(v.x(1) / result[0]); //theta
                break;
            case 3:
                result[0] = Math.sqrt(v.squared());
                result[1] = Math.acos(v.x(2) / result[0]); //theta
                result[2] = Math.atan2(v.x(2), v.x(1)); //phi
                break;
            default:
                throw new IllegalArgumentException("I'm impressed by your ability to make a "+v.getD()+" dimension vector.  Please try 1-3 next time");
        }
    }
}
