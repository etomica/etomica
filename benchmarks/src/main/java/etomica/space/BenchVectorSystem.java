package etomica.space;

import etomica.box.system.VectorSystem;
import etomica.box.system.ViewVector3D;
import etomica.space3d.Space3D;
import org.openjdk.jmh.annotations.*;
import org.openjdk.jmh.profile.*;
import org.openjdk.jmh.runner.Runner;
import org.openjdk.jmh.runner.RunnerException;
import org.openjdk.jmh.runner.options.Options;
import org.openjdk.jmh.runner.options.OptionsBuilder;

import java.util.Random;
import java.util.concurrent.TimeUnit;

@State(Scope.Benchmark)
@Measurement(iterations = 10, time = 500, timeUnit = TimeUnit.MILLISECONDS)
@Warmup(iterations = 5, time = 1)
@BenchmarkMode(Mode.AverageTime)
@OutputTimeUnit(TimeUnit.NANOSECONDS)
@Fork(value = 1)
public class BenchVectorSystem {

    VectorSystem vectorSystem;
    Vector singleVec;
    VectorSystem singleVecSys;
    double scalar;

    @Setup
    public void setUp() {
        vectorSystem = new VectorSystem(Space3D.getInstance(), 1000);
        singleVec = Vector.of(1.0, 2.0, 3.0);
        singleVecSys = new VectorSystem(Space3D.getInstance(), 1);
        scalar = new Random().nextDouble();
    }

    @Benchmark
    public Vector benchAddViewView() {
        vectorSystem.get(5).PE(vectorSystem.get(6));
        return vectorSystem.get(5);
    }

    @Benchmark
    public Vector benchAddSingleSingle() {
        singleVec.PE(singleVec);
        return singleVec;
    }

    @Benchmark
    public Vector benchAddViewSingle() {
        vectorSystem.get(5).PE(singleVec);
        return vectorSystem.get(5);
    }

    @Benchmark
    public Vector benchAddViewScalar() {
        vectorSystem.get(5).PE(scalar);
        return vectorSystem.get(5);
    }

    @Benchmark
    public Vector benchAddSingleScalar() {
        singleVec.PE(scalar);
        return singleVec;
    }

    public static void main(String[] args) throws RunnerException {
        Options opts = new OptionsBuilder()
                .include(BenchVectorSystem.class.getSimpleName())
//                .jvmArgs(
//                        "-XX:+UnlockDiagnosticVMOptions",
//                        "-XX:+TraceClassLoading",
//                        "-XX:+LogCompilation",
//                        "-XX:+PrintAssembly",
//                        "-XX:-UseCompressedOops"
//                )
                .build();

        new Runner(opts).run();
    }
}
