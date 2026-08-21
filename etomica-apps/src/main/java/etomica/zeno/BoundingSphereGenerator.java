/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.zeno;

import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.space.Vector;

public class BoundingSphereGenerator {
    public static class BoundingSphere {
        public Vector center;
        public double radius;
    }

    public static BoundingSphere getBoundingSphere(Box box) {
        BoundingSphere bs = new BoundingSphere();
        bs.center = box.getSpace().makeVector();
        bs.radius = 0;
        // center = getCenterRitter
        Vector rMin = box.getSpace().makeVector();
        Vector rMax = box.getSpace().makeVector();
        Vector r = box.getSpace().makeVector();
        for (IAtom a : box.getLeafList()) {
            double atomRadius = 0.5;
            r.E(a.getPosition());
            for (int i=0; i<rMin.getD(); i++) {
                if (r.getX(i) - atomRadius < rMin.getX(i)) rMin.setX(i, r.getX(i) - atomRadius);
                if (r.getX(i) + atomRadius > rMax.getX(i)) rMax.setX(i, r.getX(i) + atomRadius);
            }
        }
        Vector initialPoint = box.getSpace().makeVector();
        initialPoint.Ev1Pv2(rMin, rMax);
        initialPoint.TE(0.5);

        //   edgePoint1 = findFarthestPoint(initialPoint)

        Vector edgePoint1 = findFarthestPoint(box, initialPoint);
        Vector edgePoint2 = findFarthestPoint(box, edgePoint1);
        bs.center.Ev1Pv2(edgePoint1, edgePoint2);
        bs.center.TE(0.5);
        Vector edgePoint = findFarthestPoint(box, bs.center);
        bs.radius = Math.sqrt(edgePoint.Mv1Squared(bs.center));

        return bs;
    }

    public static Vector findFarthestPoint(Box box, Vector queryPoint) {
        Vector farthestPoint = box.getSpace().makeVector();
        farthestPoint.E(queryPoint);
        Vector aFarthest = box.getSpace().makeVector();
        double maxR2 = 0;
        for (IAtom a : box.getLeafList()) {
            Vector r = a.getPosition();
            if (r.Mv1Squared(queryPoint) == 0) {
                aFarthest.E(r);
                aFarthest.setX(0, aFarthest.getX(0) + 0.5);
            }
            else {
                aFarthest.Ev1Mv2(r, queryPoint);
                aFarthest.normalize();
                aFarthest.TE(0.5);
                aFarthest.PE(r);
            }
            double r2 = aFarthest.Mv1Squared(queryPoint);
            if (r2 > maxR2) {
                farthestPoint.E(aFarthest);
                maxR2 = r2;
            }
        }
        return farthestPoint;
    }
}
