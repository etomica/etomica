/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.config;

import etomica.atom.IAtomList;
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

    private static double[] parseLine(String line) {
        return Arrays.stream(line.split("\\s+"))
                .mapToDouble(Double::valueOf)
                .toArray();
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
}
