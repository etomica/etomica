/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.zeno;

import etomica.space.Vector;
import etomica.space3d.RotationTensor3D;
import etomica.space3d.Vector3D;
import etomica.util.random.IRandom;

/**
 * This reproduces the math from ZENO's class of the same name.
 */
public class BiasedSpherePointDirect {
    public static double computeCosTheta(double alpha, double R) {
        double num = -Math.pow(1.0 - alpha, 2.0) + 2.0 * (1.0 - alpha) * (1.0 + alpha*alpha) * R + 2.0 * alpha * (1.0 + alpha*alpha) * R*R;
        double den = Math.pow(1.0 - alpha + 2.0 * alpha * R, 2.0);
        double cosTheta = num / den;
        if (cosTheta > 1.0) {
            cosTheta = 1.0;
        }

        return cosTheta;
    }

    public static Vector3D generate(IRandom random, Vector sphereCenter, double sphereRadius, Vector3D distributionCenter, double alpha) {
        double R = random.nextDouble();
        double cosTheta = computeCosTheta(alpha, R);
        double sinTheta = Math.sqrt(1.0 - Math.pow(cosTheta, 2.0));
        double phi = (Math.PI * 2) * random.nextDouble();

        Vector3D point = new Vector3D(sinTheta * Math.cos(phi), sinTheta * Math.sin(phi), cosTheta);
        Vector3D recenteredDistributionCenter = new Vector3D();
        recenteredDistributionCenter.Ev1Mv2(distributionCenter, sphereCenter);
        recenteredDistributionCenter.normalize();
        Vector3D rotationAxis = new Vector3D();
        rotationAxis.E(recenteredDistributionCenter);
        rotationAxis.XE(new Vector3D(0, 0, 1));
        double sinRotationAngle = Math.sqrt(rotationAxis.squared());
        double cosRotationAngle = recenteredDistributionCenter.getX(2);
        rotationAxis.normalize();
        RotationTensor3D rotationTensor = new RotationTensor3D();
        rotationTensor.setRotationAxis(rotationAxis, -Math.atan2(sinRotationAngle, cosRotationAngle));
        rotationTensor.transform(point);
        point.TE(sphereRadius);
        point.PE(sphereCenter);
        return point;
    }
}
