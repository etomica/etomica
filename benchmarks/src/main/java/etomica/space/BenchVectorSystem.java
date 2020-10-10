package etomica.space;

import etomica.box.storage.VectorStorage;
import etomica.space3d.Space3D;
import org.openjdk.jmh.annotations.*;
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

    VectorStorage vectorStorage;
    Vector singleVec;
    VectorStorage singleVecSys;
    double scalar;

    @Setup
    public void setUp() {
        vectorStorage = new VectorStorage(Space3D.getInstance(), 1000);
        singleVec = Vector.of(1.0, 2.0, 3.0);
        singleVecSys = new VectorStorage(Space3D.getInstance(), 1);
        scalar = new Random().nextDouble();
    }

    @Benchmark
    public Vector benchAddViewView() {
        vectorStorage.get(5).PE(vectorStorage.get(6));
        return vectorStorage.get(5);
    }

    @Benchmark
    public Vector benchAddSingleSingle() {
        singleVec.PE(singleVec);
        return singleVec;
    }

    @Benchmark
    public Vector benchAddViewSingle() {
        vectorStorage.get(5).PE(singleVec);
        return vectorStorage.get(5);
    }

    @Benchmark
    public Vector benchAddViewScalar() {
        vectorStorage.get(5).PE(scalar);
        return vectorStorage.get(5);
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
