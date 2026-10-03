/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.util.voro.examples;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.InputStream;

public class ResourceHelper {

    public static InputStream getStreamForFile(String filename, Class c) {
        InputStream in = null;
        if (new File(filename).exists()) {
            // local file exists; use that
            try {
                in = new FileInputStream(filename);
            } catch (FileNotFoundException e) {
                throw new RuntimeException(e);
            }
        }
        else {
            // use file from resources
            in = c.getResourceAsStream(filename);
        }
        return in;
    }
}
