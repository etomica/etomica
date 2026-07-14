/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.zeno;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.space.Vector;
import etomica.space3d.Vector3D;
import etomica.util.random.IRandom;

public class WalkerExterior {
    protected final IRandom random;
    protected final Box box;
    protected final double[] sigmaByType;
    protected final double boundingSphereRadius;
    protected final Vector boundingSphereCenter;
    protected final double shellThickness;

    public WalkerExterior(Box box, double[] sigmaByType, IRandom random, double boundingSphereRadius, Vector boundingSphereCenter, double shellThickness) {
        this.box = box;
        this.random = random;
        this.sigmaByType = sigmaByType;
        this.boundingSphereRadius = boundingSphereRadius;
        this.boundingSphereCenter = boundingSphereCenter;
        this.shellThickness = shellThickness;
    }

    public double findNearestPoint(Vector position) {
        IAtomList atoms = this.box.getLeafList();
        double minDistance = Double.POSITIVE_INFINITY;

        for(IAtom a : atoms) {
            double r2 = a.getPosition().Mv1Squared(position);
            double sigma = this.sigmaByType[a.getType().getIndex()];
            double distance = Math.sqrt(r2) - sigma * 0.5;
            if (distance < minDistance) {
                minDistance = distance;
            }
        }

        return minDistance * minDistance;
    }

    public WalkResult walk() {
        WalkResult result = new WalkResult();
        Vector3D position = new Vector3D();
        position.setRandomSphere(this.random);
        position.TE(this.boundingSphereRadius);
        position.PE(this.boundingSphereCenter);
//        System.out.println("center "+this.boundingSphereCenter);
//        System.out.println("radius "+this.boundingSphereRadius);
//        position.E(-7.327586 , -19.809426 , 18.303763);
//        position.E( -7.327586 , -19.809426 , 18.303763);
        result.startPoint.E(position);

        while(true) {
//            System.out.println("position "+position);
            double minDistanceSqr = this.findNearestPoint(position);
//            System.out.println("minDistnaceSqr "+minDistanceSqr);
            double minDistance = Math.sqrt(minDistanceSqr);
            if (minDistance <= this.shellThickness) {
                result.endPoint.E(position);
                result.hitObject = true;
//                System.out.println("walker absorbed");
//                System.exit(0);
                return result;
            }

            ++result.numSteps;
            Vector3D step = new Vector3D();
            step.setRandomSphere(this.random);
            step.TE(minDistance);
            position.PE(step);
//            System.out.println("after step "+position);
//            if (result.numSteps == 1) {
//                position.E(-12.240174, -6.256060, 27.351211);
//                System.out.println("after step => " + position);
//            }
            double centerDistSqr = position.Mv1Squared(this.boundingSphereCenter);
//            System.out.println("centerDistSqr "+centerDistSqr);
            if (centerDistSqr > this.boundingSphereRadius * this.boundingSphereRadius) {
                double alpha = this.boundingSphereRadius / Math.sqrt(centerDistSqr);
                if (this.random.nextDouble() < 1.0 - alpha) {
                    result.hitObject = false;
//                    System.out.println("walker escapes");
//                    System.exit(0);
                    return result;
                }

                position = BiasedSpherePointDirect.generate(this.random, this.boundingSphereCenter, this.boundingSphereRadius, position, alpha);
//                System.out.println("BSPD => "+position);
            }
        }
    }

    public static class WalkResult {
        public boolean hitObject;
        public int numSteps;
        public Vector3D startPoint = new Vector3D();
        public Vector3D endPoint = new Vector3D();
    }
}

