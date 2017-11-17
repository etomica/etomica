/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.math;

import org.openjdk.jmh.annotations.*;

import java.util.concurrent.TimeUnit;


/**
 * Created by alex on 4/29/17.
 */
@State(Scope.Benchmark)
@Measurement(time = 100, timeUnit = TimeUnit.MILLISECONDS)
@BenchmarkMode(Mode.AverageTime)
@OutputTimeUnit(TimeUnit.NANOSECONDS)
public class FactorialBenchmark {

    @Param({"1", "5", "20"})
    public int n;

    @Benchmark
    public long etomicaFactorial() {
        return etomica.math.SpecialFunctions.factorial(n);
    }

    @Benchmark
    public long apacheFactorial() {
        return org.apache.commons.math3.util.CombinatoricsUtils.factorial(n);
    }
}
