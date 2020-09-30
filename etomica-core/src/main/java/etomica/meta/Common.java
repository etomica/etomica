package etomica.meta;

import io.github.classgraph.ClassGraph;
import io.github.classgraph.ScanResult;

public class Common {

    public static final ScanResult CLASSPATH_SCAN = new ClassGraph().acceptPackages("etomica").verbose().scan();

}
