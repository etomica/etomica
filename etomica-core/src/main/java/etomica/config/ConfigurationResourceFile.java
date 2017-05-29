package etomica.config;

import etomica.api.IAtomList;
import etomica.box.Box;
import etomica.space.Vector;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Arrays;

/**
 *
 */
public class ConfigurationResourceFile implements Configuration {

    private final String fileName;
    private final Class callingClass;

    /**
     *
     * @param fileName
     * @param callingClass
     */
    public ConfigurationResourceFile(String fileName, Class callingClass) {
        this.fileName = fileName;
        this.callingClass = callingClass;
    }

    @Override
    public void initializeCoordinates(Box box) {
        IAtomList leafList = box.getLeafList();
        try(BufferedReader reader =
                    new BufferedReader(new InputStreamReader(callingClass.getResourceAsStream(fileName)))) {

            for(int i = 0; i < leafList.getAtomCount(); i++) {
                Vector pos = leafList.getAtom(i).getPosition();
                pos.E(parseLine(reader.readLine()));
            }
        } catch(IOException e) {
            e.printStackTrace();
        }

    }

    private static double[] parseLine(String line) {
        return Arrays.stream(line.split("\\s+"))
                .mapToDouble(Double::valueOf)
                .toArray();
    }
}
