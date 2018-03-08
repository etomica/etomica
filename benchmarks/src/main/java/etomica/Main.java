package etomica;

import org.openjdk.jmh.profile.LinuxPerfAsmProfiler;
import org.openjdk.jmh.profile.LinuxPerfNormProfiler;
import org.openjdk.jmh.profile.StackProfiler;
import org.openjdk.jmh.runner.Runner;
import org.openjdk.jmh.runner.RunnerException;
import org.openjdk.jmh.runner.options.Options;
import org.openjdk.jmh.runner.options.OptionsBuilder;

public class Main {
    public static void main(String[] args) throws RunnerException {

        Options opts = new OptionsBuilder()
                .include("BenchSimLJMD3D")
                .jvmArgsAppend("-XX:+UnlockDiagnosticVMOptions", "-XX:+TraceClassLoading", "-XX:+LogCompilation", "-XX:+PrintAssembly")
//                .jvmArgsAppend("-XX:InlineSmallCode=2500")
                .build();

        new Runner(opts).run();
    }
}
