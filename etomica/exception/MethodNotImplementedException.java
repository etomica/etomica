package etomica.exception;

/**
 * Exception that indicates a method is defined for a class (perhaps to adhere
 * to an interface), but its implementation is not yet completed.
 * @author David Kofke
 *
 */

/* History
 * 
 * 01/25/03 (DAK) new
 */
public class MethodNotImplementedException extends RuntimeException {

	/**
	 * Constructor for MethodNotImplementedException.
	 */
	public MethodNotImplementedException() {
		this("Attempt to invoke method that exists but is not implemented (class is still under development)");
	}

	/**
	 * Constructor for MethodNotImplementedException.
	 * @param message
	 */
	public MethodNotImplementedException(String message) {
		super(message);
	}

	/**
	 * Constructor for MethodNotImplementedException.
	 * @param message
	 * @param cause
	 */
	public MethodNotImplementedException(String message, Throwable cause) {
		super(message, cause);
	}

	/**
	 * Constructor for MethodNotImplementedException.
	 * @param cause
	 */
	public MethodNotImplementedException(Throwable cause) {
		super(cause);
	}

}
