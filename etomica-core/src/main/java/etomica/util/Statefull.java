package etomica.util;

import java.io.*;

public interface Statefull {
    void saveState(Writer fw) throws IOException;
    void restoreState(BufferedReader br) throws IOException;
}
