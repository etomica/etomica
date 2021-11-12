package etomica.simulation;

import etomica.data.IData;
import etomica.modules.dcvgcmd.DCVGCMDFasterer;
import org.openjdk.jmh.annotations.*;

import java.util.concurrent.TimeUnit;

@State(Scope.Benchmark)
@Fork(1)
public class BenchSimDCVGCMD {

    private DCVGCMDFasterer sim;

    @Setup(Level.Iteration)
    public void setUp() {
        sim = new DCVGCMDFasterer();
    }

    @Benchmark
    @BenchmarkMode(Mode.Throughput)
    @OutputTimeUnit(TimeUnit.SECONDS)
    @Warmup(time = 1, iterations = 5)
    @Measurement(time = 10, timeUnit = TimeUnit.SECONDS, iterations = 5)
    public IData integratorStep() {
        sim.integratorDCV.doStep();
        return sim.accumulator1.getData();
    }
}
