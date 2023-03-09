/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.simulation;

import etomica.config.Configuration;
import etomica.config.ConfigurationResourceFile;
import etomica.potential.compute.PotentialCompute;
import etomica.tests.TestLJMC3D;
import etomica.tests.TestLJMD3D;
import etomica.tests.TestLJMD3DBrute;
import etomica.tests.TestLJMD3DNew;
import org.openjdk.jmh.annotations.*;
import org.openjdk.jmh.runner.Runner;
import org.openjdk.jmh.runner.RunnerException;
import org.openjdk.jmh.runner.options.Options;
import org.openjdk.jmh.runner.options.OptionsBuilder;

import java.util.concurrent.TimeUnit;

@State(Scope.Benchmark)
@Fork(1)
public class BenchSimLJMD3D {

    @Param({"500", "4000"})
    private int numMolecules;

    @Param({"200000"})
    private int numSteps;

    private TestLJMD3D sim;
    private PotentialCompute pm;

    private TestLJMD3DNew simNew;

    private TestLJMD3DBrute simBrute;
    private PotentialCompute pmBrute;

    @Setup(Level.Iteration)
    public void setUp() {

        Configuration config = new ConfigurationResourceFile(
                String.format("LJMC3D%d.pos", numMolecules),
                TestLJMC3D.class
        );

        {
            sim = new TestLJMD3D(numMolecules, config);
            sim.integrator.reset();
            pm = sim.integrator.getPotentialCompute();
        }

        {
            simNew = new TestLJMD3DNew(numMolecules, config);
            simNew.integrator.reset();
        }

        {
            simBrute = new TestLJMD3DBrute(numMolecules, config);
            simBrute.integrator.reset();
            pmBrute = sim.integrator.getPotentialCompute();
        }
    }

    @Benchmark
    @BenchmarkMode(Mode.Throughput)
    @OutputTimeUnit(TimeUnit.SECONDS)
    @Warmup(time = 1, iterations = 5)
    @Measurement(time = 3, iterations = 5)
    public void integratorStep() {
        sim.integrator.doStep();
//        pm.computeAll(true);
    }

    @Benchmark
    @BenchmarkMode(Mode.Throughput)
    @OutputTimeUnit(TimeUnit.SECONDS)
    @Warmup(time = 1, iterations = 5)
    @Measurement(time = 3, iterations = 5)
    public void integratorStepNew() {
        simNew.integrator.doStep();
//        pm.computeAll(true);
    }

//    @Benchmark
    @BenchmarkMode(Mode.Throughput)
    @OutputTimeUnit(TimeUnit.SECONDS)
    @Warmup(time = 1, iterations = 5)
    @Measurement(time = 3, iterations = 5)
    public void integratorStepBrute() {
        simBrute.integrator.doStep();
//        pmBrute.computeAll(true);
    }

    public static void main(String[] args) throws RunnerException {

        Options opts = new OptionsBuilder()
                .include(BenchSimLJMD3D.class.getSimpleName())
//                .jvmArgs(
//                        "-XX:+UnlockDiagnosticVMOptions",
//                        "-XX:+PrintAssembly",
//                        "-XX:PrintAssemblyOptions=intel",
//                        "-XX:CompileCommand=print,*BoxBench.bench*")
                .build();

        new Runner(opts).run();
    }
}
