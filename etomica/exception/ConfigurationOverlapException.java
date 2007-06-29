package etomica.exception;

import etomica.box.Box;

/**
 * Exception thrown when an overlap is detected in the configuration of
 * a box, subject to the status of an ignoreOverlap flag.
 */
public class ConfigurationOverlapException extends Exception {

    /**
     * Constructor for ConfigurationOverlapException
     * 
     * @param box the Box in which the overlap is detected
     */
    public ConfigurationOverlapException(Box box) {
        super("Overlap in configuration "+box);
        this.box = box;
    }

    public final Box box;
}
