package etomica.math;

import org.openjdk.jmh.annotations.Benchmark;
import org.openjdk.jmh.annotations.Param;
import org.openjdk.jmh.annotations.Scope;
import org.openjdk.jmh.annotations.State;


/**
 * Created by alex on 4/29/17.
 */
@State(Scope.Benchmark)
public class FactorialBenchmark {

    @Param({"1", "5", "20"})
    public int n;

    @Benchmark
    public long benchEtomicaFactorial() {
        return etomica.math.SpecialFunctions.factorial(n);
    }

    @Benchmark
    public long benchApacheFactorial() {
        return org.apache.commons.math3.util.CombinatoricsUtils.factorial(n);
    }
}
