package etomica.exception;

import etomica.phase.Phase;

/**
 * Exception thrown when an overlap is detected in the configuration of
 * a phase, subject to the status of an ignoreOverlap flag.
 */
public class ConfigurationOverlapException extends Exception {

    /**
     * Constructor for ConfigurationOverlapException
     * 
     * @param phase the Phase in which the overlap is detected
     */
    public ConfigurationOverlapException(Phase phase) {
        super("Overlap in configuration "+phase);
        this.phase = phase;
    }

    public final Phase phase;
}
