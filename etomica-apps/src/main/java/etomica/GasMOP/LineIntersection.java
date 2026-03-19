package etomica.GasMOP;

import etomica.space.Vector;
import etomica.space3d.Vector3D;

public class LineIntersection {
    public Vector findIntersection(Line line1, Line line2, double tolerance) {
        // Parametric equations for the two lines:
        // L1(t) = P1 + t*V1
        // L2(s) = P2 + s*V2
        // We set L1(t) = L2(s) and solve for t and s.

        Vector p1 = line1.point;
        Vector v1 = line1.direction;
        Vector p2 = line2.point;
        Vector v2 = line2.direction;

        // The system of equations to solve is:
        // P1.x + t*V1.x = P2.x + s*V2.x  =>  t*V1.x - s*V2.x = P2.x - P1.x (Eq 1)
        // P1.y + t*V1.y = P2.y + s*V2.y  =>  t*V1.y - s*V2.y = P2.y - P1.y (Eq 2)
        // P1.z + t*V1.z = P2.z + s*V2.z  =>  t*V1.z - s*V2.z = P2.z - P1.z (Eq 3)

        // We solve the first two equations for t and s using Cramer's Rule or substitution.
        // We'll choose the two equations with the largest determinant to avoid division by zero.
        /*double detXY = v1.x * (-v2.y) - (-v2.x) * v1.y;
        double detYZ = v1.y * (-v2.z) - (-v2.y) * v1.z;
        double detXZ = v1.x * (-v2.z) - (-v2.x) * v1.z;*/

        double detXY = v1.getX(0) * (-v2.getX(1)) - (-v2.getX(0)) * v1.getX(1);
        double detYZ = v1.getX(1) * (-v2.getX(2)) - (-v2.getX(1)) * v1.getX(2);
        double detXZ = v1.getX(0) * (-v2.getX(2)) - (-v2.getX(0)) * v1.getX(2);

        double det = 0;
        int primaryEq1 = -1, primaryEq2 = -1;

        if (Math.abs(detXY) > Math.abs(detYZ) && Math.abs(detXY) > Math.abs(detXZ)) {
            det = detXY;
            primaryEq1 = 1;
            primaryEq2 = 2;
        } else if (Math.abs(detYZ) > Math.abs(detXZ)) {
            det = detYZ;
            primaryEq1 = 2;
            primaryEq2 = 3;
        } else {
            det = detXZ;
            primaryEq1 = 1;
            primaryEq2 = 3;
        }

        // If the determinant is close to zero, the lines are parallel or coincide.
        if (Math.abs(det) < 1e-9) {
            // If the lines are parallel, check if they are the same line
            Vector3D diff =new Vector3D();
            diff.Ev1Mv2(p2, p1);
            if (Math.abs(v1.getX(0) * diff.getX(1) - diff.getX(0) * v1.getX(1)) < 1e-9 &&
                    Math.abs(v1.getX(0) * diff.getX(2) - diff.getX(0) * v1.getX(2)) < 1e-9) {
                // Lines are collinear (same line) - infinite intersections
                // For this example, we'll return null to indicate no unique point.
                return null;
            }
            return null; // Parallel and non-intersecting
        }

        double dx = p2.getX(0) - p1.getX(0);
        double dy = p2.getX(1) - p1.getX(1);
        double dz = p2.getX(2) - p1.getX(2);

        double t, s;

        if (primaryEq1 == 1 && primaryEq2 == 2) {
            t = (dx * (-v2.getX(1)) - dy * (-v2.getX(0))) / det;
            s = (v1.getX(0) * dy - v1.getX(1) * dx) / det;
        } else if (primaryEq1 == 2 && primaryEq2 == 3) {
            t = (dy * (-v2.getX(2)) - dz * (-v2.getX(1))) / det;
            s = (v1.getX(1) * dz - v1.getX(2) * dy) / det;
        } else { // primaryEq1 == 1 && primaryEq2 == 3
            t = (dx * (-v2.getX(2)) - dz * (-v2.getX(0))) / det;
            s = (v1.getX(0) * dz - v1.getX(2) * dx) / det;
        }

        // Check for consistency with the third, unused equation.
        double checkX = p1.getX(0) + t * v1.getX(0);
        double checkY = p1.getX(1) + t * v1.getX(1);
        double checkZ = p1.getX(2) + t * v1.getX(2);

        double checkX2 = p2.getX(0) + s * v2.getX(0);
        double checkY2 = p2.getX(1) + s * v2.getX(1);
        double checkZ2 = p2.getX(2) + s * v2.getX(2);
        Vector3D prd = new Vector3D();
        prd.Ea1Tv1(t, v1);
        prd.PE(p1);
        if (Math.abs(checkX - checkX2) < tolerance &&
                Math.abs(checkY - checkY2) < tolerance &&
                Math.abs(checkZ - checkZ2) < tolerance) {

            // The intersection point exists. Return the point using either t or s.
            return prd;
        }
        return prd;

        //return null; // Lines are skew
    }

    public Vector findPointAtDistance(Vector pointA, Vector pointB, double distance){
        Vector directionVector = new Vector3D();
        directionVector.Ev1Mv2(pointA, pointB);
        Vector directionVectorNorm = new Vector3D();
        directionVectorNorm.E(directionVector);
        directionVectorNorm.normalize();
        Vector finalVect = new Vector3D();
        finalVect.XE(directionVectorNorm);
        finalVect.TE(distance);
        return finalVect;
    }
}
