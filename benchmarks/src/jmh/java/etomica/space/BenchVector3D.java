/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.space;

import etomica.space3d.Vector3D;
import org.openjdk.jmh.annotations.*;
import org.openjdk.jmh.runner.Runner;
import org.openjdk.jmh.runner.RunnerException;
import org.openjdk.jmh.runner.options.Options;
import org.openjdk.jmh.runner.options.OptionsBuilder;

import java.util.concurrent.TimeUnit;

@State(Scope.Benchmark)
@Measurement(iterations = 10, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Warmup(iterations = 5, time = 1)
@BenchmarkMode(Mode.AverageTime)
@OutputTimeUnit(TimeUnit.NANOSECONDS)
@Fork(value = 1)
public class BenchVector3D {

    private double[] v1arr = new double[] {-30.2, 70, 20.1};
    private double[] v2arr = new double[] {50, 50, 50};

    private Vector3D v1;
    private Vector3D v2;

    @Setup
    public void setUp() {
        v1 = new Vector3D(v1arr);
        v2 = new Vector3D(v2arr);

    }

    @Benchmark
    public double baseline() {
        return v1.getX(0) + v1.getX(1) + v1.getX(2);
    }

    @Benchmark
    public double mod() {
        v1.mod(v2);

        return v1.getX(0) + v1.getX(2) + v1.getX(1);
    }

    public static void main(String[] args) throws RunnerException {
        Options opts = new OptionsBuilder()
                .include(BenchVector3D.class.getSimpleName())
                .build();

        new Runner(opts).run();

    }

}
