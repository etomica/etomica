package etomica.space;

import etomica.space3d.Vector3D;
import org.openjdk.jmh.annotations.*;
import org.openjdk.jmh.runner.Runner;
import org.openjdk.jmh.runner.RunnerException;
import org.openjdk.jmh.runner.options.Options;
import org.openjdk.jmh.runner.options.OptionsBuilder;

import java.util.Random;
import java.util.concurrent.TimeUnit;

@State(Scope.Benchmark)
@Measurement(iterations = 20, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Warmup(iterations = 10, time = 1)
@BenchmarkMode(Mode.SingleShotTime)
@OutputTimeUnit(TimeUnit.MILLISECONDS)
@Fork(value = 1)
public class BenchVector {

    private final Vector3D[] drs = new Vector3D[10_000_000];
    private final Vector3D dim = new Vector3D(36.0, 36.0, 36.0);
    private final Vector3D dimHalf = new Vector3D();

    {
        for (int i = 0; i < drs.length; i++) {
            drs[i] = new Vector3D();
        }
    }

    @Setup
    public void setUp() {
        dimHalf.Ea1Tv1(0.5, dim);
        Random r = new Random();
        for (Vector dr : drs) {
            dr.E(r.doubles(3, -70.0, 70.0).toArray());
        }
    }

    @Benchmark
    public Vector[] nearestImageOriginal() {
        for (Vector3D dr : drs) {
            dr.PE(dimHalf);
            dr.mod(dim);
            dr.ME(dimHalf);
        }
        return drs;
    }

    @Benchmark
    public Vector[] nearestImageOpt() {
        for (Vector3D dr : drs) {
            dr.nearestImage(dim);
        }
        return drs;
    }


    public static void main(String[] args) throws RunnerException {
        Options opts = new OptionsBuilder()
                .include(BenchVector.class.getSimpleName())
                .exclude(BenchVector3D.class.getSimpleName())
                .build();

        new Runner(opts).run();

    }
}
