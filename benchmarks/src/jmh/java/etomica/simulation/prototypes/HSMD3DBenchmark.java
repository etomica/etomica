/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.simulation.prototypes;

import etomica.simulation.Simulation;
import org.openjdk.jmh.annotations.Benchmark;
import org.openjdk.jmh.annotations.Scope;
import org.openjdk.jmh.annotations.Setup;
import org.openjdk.jmh.annotations.State;
import org.openjdk.jmh.profile.GCProfiler;
import org.openjdk.jmh.runner.Runner;
import org.openjdk.jmh.runner.RunnerException;
import org.openjdk.jmh.runner.options.Options;
import org.openjdk.jmh.runner.options.OptionsBuilder;

/**
 * Created by alex on 4/29/17.
 */
@State(Scope.Benchmark)
public class HSMD3DBenchmark {

    private Simulation sim;

    @Setup
    public void setUp() {
        sim = new HSMD3D();
    }

    @Benchmark
    public void hsmd3d() {
        sim.getIntegrator().doStep();
    }

    public static void main(String[] args) throws RunnerException {
        Options opts = new OptionsBuilder()
                .include(HSMD3DBenchmark.class.getSimpleName())
                .forks(0)
                .addProfiler(GCProfiler.class)
                .build();

        new Runner(opts).run();
    }
}
