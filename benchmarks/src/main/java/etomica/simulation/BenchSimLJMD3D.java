/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.simulation;

import etomica.config.Configuration;
import etomica.config.Configurations;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.IteratorDirective;
import etomica.potential.PotentialCalculationForceSum;
import etomica.simulation.prototypes.LJMD3D;
import etomica.simulation.prototypes.LJMD3DNbr;
import etomica.tests.TestLJMC3DSlowerer;
import etomica.tests.TestLJMD3D;
import etomica.tests.TestLJMD3DSlowerer;
import org.openjdk.jmh.annotations.*;
import org.openjdk.jmh.profile.LinuxPerfAsmProfiler;
import org.openjdk.jmh.profile.LinuxPerfNormProfiler;
import org.openjdk.jmh.profile.Profiler;
import org.openjdk.jmh.runner.Runner;
import org.openjdk.jmh.runner.RunnerException;
import org.openjdk.jmh.runner.options.Options;
import org.openjdk.jmh.runner.options.OptionsBuilder;

import java.util.concurrent.TimeUnit;

@State(Scope.Benchmark)
@Fork(1)
public class BenchSimLJMD3D {

    private TestLJMD3DSlowerer sim;
    private TestLJMD3D simFast;

    @Param({"500", "4000"})
    private int numAtoms;

    @Param({"200000"})
    private int numSteps;

    @Setup(Level.Iteration)
    public void setUp() {
        TestLJMD3D.SimParams params = new TestLJMD3D.SimParams();
        params.numAtoms = numAtoms;
        params.numSteps = numSteps;
        Configuration config = Configurations.fromResourceFile(String.format("LJMC3D%d.pos", params.numAtoms), TestLJMC3DSlowerer.class);

        sim = new TestLJMD3DSlowerer(params.numAtoms, config);
        simFast = new TestLJMD3D(params.numAtoms, config);
        sim.integrator.reset();
        simFast.integrator.reset();

    }

    @Benchmark
    @BenchmarkMode(Mode.Throughput)
    @OutputTimeUnit(TimeUnit.SECONDS)
    @Warmup(time = 1, iterations = 5)
    @Measurement(time = 3, iterations = 5)
    public long integratorStep() {
        sim.integrator.doStep();
//        pm.calculate(sim.box, pc, false);
        return sim.integrator.getStepCount();
    }

    @Benchmark
    @BenchmarkMode(Mode.Throughput)
    @OutputTimeUnit(TimeUnit.SECONDS)
    @Warmup(time = 1, iterations = 5)
    @Measurement(time = 3, iterations = 5)
    public long integratorStepFast() {
        simFast.integrator.doStep();
//        pm.calculate(sim.box, pc, false);
        return simFast.integrator.getStepCount();
    }


    public static void main(String[] args) throws RunnerException {

        Options opts = new OptionsBuilder()
                .include(BenchSimLJMD3D.class.getSimpleName())
//                .jvmArgs(
//                        "-XX:+UnlockDiagnosticVMOptions",
//                        "-XX:+PrintAssembly",
//                        "-XX:PrintAssemblyOptions=intel",
//                        "-XX:CompileCommand=print,*BoxBench.bench*")
                .addProfiler(LinuxPerfAsmProfiler.class, "frequency=5000;intelSyntax=true")
//                .addProfiler(LinuxPerfNormProfiler.class)
                .build();

        new Runner(opts).run();
    }
}
