/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.simulation;

import etomica.config.Configuration;
import etomica.config.ConfigurationResourceFile;
import etomica.potential.PotentialCalculationForceSum;
import etomica.potential.PotentialMaster;
import etomica.potential.compute.PotentialCompute;
import etomica.tests.*;
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

    private TestLJMD3DBrute simBrute;
    private PotentialCompute pmBrute;

    private TestLJMD3DSlowerer simSlow;
    private PotentialCalculationForceSum pcSlow;
    private PotentialMaster pmSlow;

    private TestLJMD3DSlowBrute simSlowBrute;
    private PotentialCalculationForceSum pcSlowBrute;
    private PotentialMaster pmSlowBrute;

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
            simBrute = new TestLJMD3DBrute(numMolecules, config);
            simBrute.integrator.reset();
            pmBrute = sim.integrator.getPotentialCompute();
        }

        {
            simSlow = new TestLJMD3DSlowerer(numMolecules, config);
            simSlow.integrator.reset();
            pcSlow = simSlow.integrator.getForceSum();
            pmSlow = simSlow.integrator.getPotentialMaster();
        }


        {
            simSlowBrute = new TestLJMD3DSlowBrute(numMolecules, config);
            simSlowBrute.integrator.reset();
            pcSlowBrute = simSlowBrute.integrator.getForceSum();
            pmSlowBrute = simSlowBrute.integrator.getPotentialMaster();
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
    public void integratorStepBrute() {
        simBrute.integrator.doStep();
//        pmBrute.computeAll(true);
    }

    @Benchmark
    @BenchmarkMode(Mode.Throughput)
    @OutputTimeUnit(TimeUnit.SECONDS)
    @Warmup(time = 1, iterations = 5)
    @Measurement(time = 3, iterations = 5)
    public void integratorStepSlowerer() {
        simSlow.integrator.doStep();
//        pmSlow.calculate(simSlow.box, null, pcSlow);
    }

    @Benchmark
    @BenchmarkMode(Mode.Throughput)
    @OutputTimeUnit(TimeUnit.SECONDS)
    @Warmup(time = 1, iterations = 5)
    @Measurement(time = 3, iterations = 5)
    public void integratorStepSlowBrute() {
        simSlowBrute.integrator.doStep();
//        pmSlowBrute.calculate(simSlowBrute.box, null, pcSlowBrute);
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
