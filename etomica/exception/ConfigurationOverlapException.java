package etomica.exception;

import etomica.phase.Phase;

public class ConfigurationOverlapException extends Exception {

    /**
     * Constructor for MethodNotImplementedException.
     */
    public ConfigurationOverlapException(Phase phase) {
	    super("Overlap in configuration "+phase);
	    this.phase = phase;
    }

    public Phase phase;
}
