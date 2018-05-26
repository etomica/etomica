/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.simulation;

import etomica.nbr.list.PotentialMasterList;
import etomica.potential.IteratorDirective;
import etomica.potential.PotentialCalculationForceSum;
import etomica.simulation.prototypes.LJMD3D;
import etomica.simulation.prototypes.LJMD3DNbr;
import org.openjdk.jmh.annotations.*;
import org.openjdk.jmh.runner.Runner;
import org.openjdk.jmh.runner.RunnerException;
import org.openjdk.jmh.runner.options.Options;
import org.openjdk.jmh.runner.options.OptionsBuilder;

import java.util.concurrent.TimeUnit;

@State(Scope.Benchmark)
@Fork(1)
public class BenchSimLJMD3D {

    private LJMD3DNbr sim;
    private PotentialCalculationForceSum pc;
    private IteratorDirective id;
    private PotentialMasterList pm;

    @Setup(Level.Iteration)
    public void setUp() {

        sim = new LJMD3DNbr();
        sim.integrator.reset();
        id = new IteratorDirective(IteratorDirective.Direction.UP);
        pc = sim.integrator.getForceSum();
        pm = (PotentialMasterList) sim.integrator.getPotentialMaster();

    }

    @Benchmark
    @BenchmarkMode(Mode.Throughput)
    @OutputTimeUnit(TimeUnit.SECONDS)
    @Warmup(time = 1, iterations = 5)
    @Measurement(time = 3, iterations = 5)
    public long integratorStep() {
//        sim.integrator.doStep();
        pm.calculate(sim.box, pc, false);
        return sim.integrator.getStepCount();
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
