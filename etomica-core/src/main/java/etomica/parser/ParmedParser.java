package etomica.parser;

import java.io.IOException;
import java.io.InputStreamReader;

/**
 * Created by alex on 3/16/17.
 */
public class ParmedParser {
    public static void main(String[] args) {
        try {
            ProcessBuilder pb = new ProcessBuilder("venv/bin/parmed_json", "test.top", "test.gro");
            System.out.println(pb.environment());
            pb.inheritIO();

            Process p = pb.start();
            p.waitFor();

        } catch (IOException | InterruptedException e) {
            e.printStackTrace();
        }
    }
}
