package etomica.math;

import org.openjdk.jmh.annotations.*;
import org.openjdk.jmh.runner.Runner;
import org.openjdk.jmh.runner.RunnerException;
import org.openjdk.jmh.runner.options.Options;
import org.openjdk.jmh.runner.options.OptionsBuilder;

import java.util.Random;
import java.util.concurrent.TimeUnit;

@State(Scope.Benchmark)
@Warmup(iterations = 3, time = 3)
@Measurement(iterations = 5, time = 3)
@BenchmarkMode(Mode.Throughput)
@Fork(1)
public class BenchComplex {

    private Complex[] eIn1;
    private Complex[] eIn2;
    private Complex[] eOut;
    private double[] eDIn1;
    private double[] eDIn2;
    private double[] eDOut;
    private org.apache.commons.math3.complex.Complex[] aIn1;
    private org.apache.commons.math3.complex.Complex[] aIn2;
    private org.apache.commons.math3.complex.Complex[] aOut;
    private double[] aDIn1;
    private double[] aDIn2;
    private double[] aDOut;

    @Setup(Level.Iteration)
    public void setUp() {
        int n = 100000;

        eIn1 = new Complex[n];
        eIn2 = new Complex[n];
        eOut = new Complex[n];
        eDIn1 = new double[n * 2];
        eDIn2 = new double[n * 2];
        eDOut = new double[n * 2];


        aIn1 = new org.apache.commons.math3.complex.Complex[n];
        aIn2 = new org.apache.commons.math3.complex.Complex[n];
        aOut = new org.apache.commons.math3.complex.Complex[n];
        aDIn1 = new double[n * 2];
        aDIn2 = new double[n * 2];
        aDOut = new double[n * 2];

        Random r = new Random();
        for (int i = 0; i < n; i++) {
            eIn1[i] = new Complex(r.nextDouble(), r.nextDouble());
            eIn2[i] = new Complex(r.nextDouble(), r.nextDouble());
            eOut[i] = Complex.ZERO;

            aIn1[i] = new org.apache.commons.math3.complex.Complex(r.nextDouble(), r.nextDouble());
            aIn2[i] = new org.apache.commons.math3.complex.Complex(r.nextDouble(), r.nextDouble());
            aOut[i] = org.apache.commons.math3.complex.Complex.ZERO;
        }

        for (int i = 0; i < eDIn1.length; i++) {
            eDIn1[i] = r.nextDouble();
            eDIn2[i] = r.nextDouble();
            eDOut[i] = r.nextDouble();

            aDIn1[i] = r.nextDouble();
            aDIn2[i] = r.nextDouble();
            aDOut[i] = r.nextDouble();
        }

    }

    @Benchmark
    public Complex[] benchEtomica() {
        for (int i = 0; i < eIn1.length; i++) {
            eOut[i] = eIn1[i].times(eIn2[i]);
        }
        return eOut;
    }

    @Benchmark
    public double[] benchEtomicaDoubles() {
        for (int i = 0; i < eIn1.length; i++) {
            Complex.fromArray(eDIn1, i).times(Complex.fromArray(eDIn2, i)).intoArray(eDOut, i);
        }
        return eDOut;
    }

    @Benchmark
    public org.apache.commons.math3.complex.Complex[] benchApache() {
        for (int i = 0; i < aIn1.length; i++) {
            aOut[i] = aIn1[i].multiply(aIn2[i]);
        }
        return aOut;
    }

    public static void main(String[] args) throws RunnerException {

        Options opts = new OptionsBuilder()
                .include("BenchComplex.benchEtomica")
                .jvmArgsAppend("-XX:-UseSuperWord", "-XX:+UnlockDiagnosticVMOptions", "-XX:CompileCommand=print,etomica/math/BenchComplex.benchEtomica")
                .build();

        new Runner(opts).run();
    }
}
