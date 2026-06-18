package etomica.util;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.Writer;

/**
 * Interface for a class that can be saved or restored.  The saveState method can be called on many objects.  In this
 * case, restoreState should be called in the same order when restoring.
 *
 * The nature of what is written (text vs. binary, formatting, etc) is entirely up to the implementation of this class.
 */
public interface Statefull {

    /**
     * Writes the state of the object via the given Writer.
     * @param fw
     * @throws IOException
     */
    void saveState(Writer fw) throws IOException;

    /**
     * Restores the state of the object from the given BufferedReader.
     * @param br
     * @throws IOException
     */
    void restoreState(BufferedReader br) throws IOException;
}
