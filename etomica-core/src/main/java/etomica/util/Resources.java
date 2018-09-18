package etomica.util;

import java.io.BufferedReader;
import java.io.InputStreamReader;

public final class Resources {
    private Resources() {}

    public static BufferedReader openResourceFile(String fileName, Class<?> cls) {
        return new BufferedReader(new InputStreamReader(cls.getResourceAsStream(fileName)));
    }
}
