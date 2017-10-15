/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.parser;

import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;

import java.io.File;
import java.io.IOException;


/**
 * Class that generates {@link ParmedStructure} objects by invoking
 * the <a href="https://github.com/ParmEd/ParmEd">ParmEd</a> library and serializing its
 * Structure object using JSON.
 *
 * @see ParmedStructure
 */
public class ParmedParser {
    private static final ObjectMapper mapper = new ObjectMapper();

    /**
     * Parses the given <a href="http://www.gromacs.org/">Gromacs</a> .top and .gro files
     * using the <a href="https://github.com/ParmEd/ParmEd">ParmEd</a> python library.
     * @param topFile File object containing the path to a Gromacs .top file
     * @param groFile File object containing the path to a Gromacs .gro file
     * @return a {@link ParmedStructure} for extracting Etomica simulation components from the ParmEd {@code Structure} object
     */
    public static ParmedStructure parseGromacs(File topFile, File groFile) throws IOException {
        JsonNode root = execParmedPython(topFile, groFile);
        return new ParmedStructure(root);
    }

    /**
     * Convenience method to parse gromacs files located in the {@code etomica-core/src/main/resources} directory.
     * @param topFileName the plain name of the .top file, e.g. {@code "test.top"}
     * @param groFileName the plain name of the .gro file, e.g. {@code "test.gro"}
     * @return a {@link ParmedStructure} containing the parsed data.
     * @throws IOException if the files do not exist.
     */
    public static ParmedStructure parseGromacsResourceFiles(String topFileName, String groFileName) throws IOException {
        return parseGromacs(getResourceFile(topFileName), getResourceFile(groFileName));
    }

    //TODO: method to accept file contents as strings and create temp files for parsing

    private static File getResourceFile(String filename) {
        ClassLoader classLoader = ParmedParser.class.getClassLoader();
        return new File(classLoader.getResource(filename).getFile());
    }

    /**
     * Runs the python script on the given files
     * @param topFile File object containing the .top file path
     * @param groFile File object containing the .gro file path
     * @return Jackson JsonNode representing the root of the json tree
     * @throws IOException if the given files do not exist
     */
    private static JsonNode execParmedPython(File topFile, File groFile) throws IOException {
        ClassLoader classLoader = ParmedParser.class.getClassLoader();
        String python_script_path = classLoader.getResource("virtualenv/bin/parmed_json").getFile();


        ProcessBuilder pb = new ProcessBuilder(
                python_script_path,
                topFile.getCanonicalPath(),
                groFile.getCanonicalPath()
        );

        Process proc = pb.start();

        JsonNode root = mapper.readTree(proc.getInputStream());

        try {
            proc.waitFor();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }

        return root;
    }
}
