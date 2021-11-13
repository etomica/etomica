/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.simulation.prototypes;

import etomica.data.AccumulatorAverageFixed;
import etomica.data.IData;
import etomica.data.types.DataDouble;
import etomica.space.Vector;
import etomica.space3d.Vector3D;
import etomica.units.Kelvin;
import etomica.units.dimensions.Null;
import etomica.util.random.IRandom;
import etomica.util.random.RandomMersenneTwister;
import etomica.util.random.RandomNumberGeneratorUnix;

public class Foo {
    public static void main(String[] args) {
        int n = 72*8;
        long steps = 100000;
        Vector v = new Vector3D();
        double mass = 120;
        IRandom random = new RandomMersenneTwister(RandomNumberGeneratorUnix.getRandSeedArray());
        AccumulatorAverageFixed acc = new AccumulatorAverageFixed(1);
        double temperatureK = 500;
        double temperature = Kelvin.UNIT.toSim(temperatureK);
        DataDouble data = new DataDouble();
        acc.putDataInfo(new DataDouble.DataInfoDouble("foo", Null.DIMENSION));
        for (int step=0; step<steps; step++) {
            double KE = 0;
            for (int a=0; a<n; a++) {
                int D = v.getD();
                for (int i = 0; i < D; i++) {
                    v.setX(i, random.nextGaussian());
                }
                v.TE(Math.sqrt(temperature / mass));
                KE += 0.5 * mass * v.squared();
            }
            double T = KE*2.0/3.0/n;
            data.x = T;
            acc.addData(data);
        }
        IData out = acc.getData();
        double avgT = Kelvin.UNIT.fromSim(out.getValue(acc.AVERAGE.index));
        double errT = Kelvin.UNIT.fromSim(out.getValue(acc.ERROR.index));
        double sdevT = Kelvin.UNIT.fromSim(out.getValue(acc.STANDARD_DEVIATION.index));
        System.out.println("avg: "+avgT+"   err: "+errT+"   sdev: "+sdevT);
        System.out.println("relative: avg: "+avgT/temperatureK+"   err: "+errT/temperatureK+"   sdev: "+sdevT/temperatureK);
    }
}
