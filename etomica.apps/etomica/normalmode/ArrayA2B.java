/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectOutputStream;


/**
 * Reads a file containing n mxm matrices as ASCII and serializes them to
 * another file.  The resulting file should be much smaller than the original
 * ASCII file.
 * 
 * This class relies on ArrayReader2D to read in the original ASCII file.
 * 
 * @author Andrew Schultz
 */
public class ArrayA2B {

    public static void main(String[] args) {
        if (args.length < 1) {
            throw new RuntimeException("usage: A2B filename");
        }
        double[][][] x = ArrayReader2D.getFromFile(args[0]);
        try {
            FileOutputStream fos = null;
            ObjectOutputStream out = null;
            fos = new FileOutputStream(args[0]+".bin");
            out = new ObjectOutputStream(fos); 
            out.writeObject(x);
            out.close();
            fos.close();
        }
        catch (IOException e) {
            throw new RuntimeException(e);
        }
    }
}
