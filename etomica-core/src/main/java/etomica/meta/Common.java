package etomica.meta;

import io.github.lukehutch.fastclasspathscanner.FastClasspathScanner;
import io.github.lukehutch.fastclasspathscanner.scanner.ScanResult;

public class Common {

    public static final ScanResult CLASSPATH_SCAN = new FastClasspathScanner("etomica").verbose().scan();

}
