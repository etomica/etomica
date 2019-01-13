/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.simulation;

import etomica.config.Configuration;
import etomica.config.ConfigurationResourceFile;
import etomica.potential.PotentialCalculationForceSum;
import etomica.potential.PotentialMaster;
import etomica.potential.PotentialMasterFasterer;
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

    @Param({"512"})
    private int totalAtoms;

    @Param({"2", "4", "8"})
    private int moleculeSize;

    private TestLJMDDimerFast simFastBrute;


    @Setup(Level.Iteration)
    public void setUp() {
        simFastBrute = new TestLJMDDimerFast(moleculeSize, totalAtoms);

    }

    @Benchmark
    @BenchmarkMode(Mode.Throughput)
    @OutputTimeUnit(TimeUnit.SECONDS)
    @Warmup(time = 1, iterations = 5)
    @Measurement(time = 3, iterations = 5)
    public void integratorStepFastBrute() {
        simFastBrute.integrator.doStep();
    }


    public static void main(String[] args) throws RunnerException {

        Options opts = new OptionsBuilder()
                .include(BenchSimLJMDDimer.class.getSimpleName())
                .addProfiler(StackProfiler.class)
//                .jvmArgs(
//                        "-XX:+UnlockDiagnosticVMOptions",
//                        "-XX:+PrintAssembly",
//                        "-XX:PrintAssemblyOptions=intel",
//                        "-XX:CompileCommand=print,*BoxBench.bench*")
                .build();

        new Runner(opts).run();
    }
}
