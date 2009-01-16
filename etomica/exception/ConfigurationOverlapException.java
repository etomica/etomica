package etomica.exception;

import etomica.api.IBox;

/**
 * Exception thrown when an overlap is detected in the configuration of
 * a box, subject to the status of an ignoreOverlap flag.
 */
public class ConfigurationOverlapException extends RuntimeException {

    /**
     * Constructor for ConfigurationOverlapException
     * 
     * @param box the Box in which the overlap is detected
     */
    public ConfigurationOverlapException(IBox box) {
        super("Overlap in configuration "+box);
        this.box = box;
    }

    private static final long serialVersionUID = 1L;
    public final IBox box;
}
