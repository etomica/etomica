package etomica.parser;

import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;
import etomica.box.Box;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;

/**
 * Created by alex on 3/16/17.
 */
public class ParmedParser {
    private static final ObjectMapper mapper = new ObjectMapper();

    public static void main(String[] args) {
        try {
            File topFile = getResourceFile("test.top");
            File groFile = getResourceFile("test.gro");


            ProcessBuilder pb = new ProcessBuilder(
                    "venv/bin/parmed_json",
                    topFile.getCanonicalPath(),
                    groFile.getCanonicalPath()
            );


            Process p = pb.start();

            JsonNode tree = mapper.readTree(p.getInputStream());

            p.waitFor();

            System.out.println(tree);

        } catch (IOException | InterruptedException e) {
            e.printStackTrace();
        }
    }

    private static File getResourceFile(String filename) {
        ClassLoader classLoader = ParmedParser.class.getClassLoader();
        return new File(classLoader.getResource(filename).getFile());
    }
}
