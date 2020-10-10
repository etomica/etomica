package etomica;

import etomica.simulation.BenchSimDCVGCMD;
import etomica.simulation.BenchSimSWChain;
import org.openjdk.jmh.profile.GCProfiler;
import org.openjdk.jmh.profile.StackProfiler;
import org.openjdk.jmh.runner.Runner;
import org.openjdk.jmh.runner.RunnerException;
import org.openjdk.jmh.runner.options.Options;
import org.openjdk.jmh.runner.options.OptionsBuilder;

public class Main {
    public static void main(String[] args) throws RunnerException {

        Options opts = new OptionsBuilder()
                .include("BenchBox.benchNewBox")
                .jvmArgs(
                        "-XX:+UnlockDiagnosticVMOptions",
                        "-XX:+TraceClassLoading",
                        "-XX:+LogCompilation",
                        "-XX:+PrintAssembly",
                        "-XX:-UseCompressedOops"
                )
//                .jvmArgs(
//                        "-XX:+UnlockDiagnosticVMOptions",
////                        "-XX:+PrintAssembly",
//                        "-XX:PrintAssemblyOptions=intel",
//                        "-XX:CompileCommand=print,etomica/BenchBox.bench*"
//                )
                .build();

        new Runner(opts).run();
    }
}
