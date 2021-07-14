/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.simulation;

import etomica.tests.*;
import org.openjdk.jmh.annotations.*;
import org.openjdk.jmh.profile.StackProfiler;
import org.openjdk.jmh.runner.Runner;
import org.openjdk.jmh.runner.RunnerException;
import org.openjdk.jmh.runner.options.Options;
import org.openjdk.jmh.runner.options.OptionsBuilder;

import java.util.concurrent.TimeUnit;

@State(Scope.Benchmark)
@Fork(1)
public class BenchSimLJMDDimer {

    @Param({"true", "false"})
    private boolean doNbr;

    @Param({"512"})
    private int totalAtoms;

    @Param({"2", "4", "16"})
    private int moleculeSize;

    private TestLJMDDimerFast simFast;
    private TestLJMDDimer simSlow;


    @Setup(Level.Iteration)
    public void setUp() {
        simFast = new TestLJMDDimerFast(moleculeSize, totalAtoms, doNbr);
        simSlow = new TestLJMDDimer(moleculeSize, totalAtoms, doNbr);

    }

    @Benchmark
    @BenchmarkMode(Mode.Throughput)
    @OutputTimeUnit(TimeUnit.SECONDS)
    @Warmup(time = 1, iterations = 5)
    @Measurement(time = 3, iterations = 5)
    public void integratorStepFast() {
        simFast.integrator.doStep();
    }

    @Benchmark
    @BenchmarkMode(Mode.Throughput)
    @OutputTimeUnit(TimeUnit.SECONDS)
    @Warmup(time = 1, iterations = 5)
    @Measurement(time = 3, iterations = 5)
    public void integratorStepSlow() {
        simSlow.integrator.doStep();
    }


    public static void main(String[] args) throws RunnerException {

        Options opts = new OptionsBuilder()
                .include(BenchSimLJMDDimer.class.getSimpleName())
//                .addProfiler(StackProfiler.class)
//                .jvmArgs(
//                        "-XX:+UnlockDiagnosticVMOptions",
//                        "-XX:+PrintAssembly",
//                        "-XX:PrintAssemblyOptions=intel",
//                        "-XX:CompileCommand=print,*BoxBench.bench*")
                .build();

        new Runner(opts).run();
    }
}
