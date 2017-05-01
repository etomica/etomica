package etomica.math;

import org.openjdk.jmh.annotations.*;

import java.util.concurrent.TimeUnit;

/**
 * Created by alex on 4/30/17.
 */
@State(Scope.Benchmark)
@Measurement(time = 100, timeUnit = TimeUnit.MILLISECONDS)
@BenchmarkMode(Mode.AverageTime)
@OutputTimeUnit(TimeUnit.NANOSECONDS)
public class LnFactorialBenchmark {

    @Param({"0", "20", "100", "10000"})
    public int n;

    @Benchmark
    public double etomicaLnFactorial() {
        return etomica.math.SpecialFunctions.lnFactorial(n);
    }

    @Benchmark
    public double apacheLnFactorial() {
        return org.apache.commons.math3.util.CombinatoricsUtils.factorialLog(n);
    }
}
