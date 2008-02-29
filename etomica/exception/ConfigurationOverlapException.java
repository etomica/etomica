package etomica.exception;

import etomica.api.IBox;

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
    public ConfigurationOverlapException(IBox box) {
        super("Overlap in configuration "+box);
        this.box = box;
    }

    public final IBox box;
}
