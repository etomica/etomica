package etomica.exception;

/**
 * @author David Kofke
 *
 * Exception that indicates a method is defined for a class (perhaps to adhere
 * to an interface), but its implementation is not yet completed.
 */

/* History
 * 
 * 01/25/03 (DAK) new
 */
public class AbnormalCompletionException extends RuntimeException {

	/**
	 * Constructor for AbnormalCompletionException
	 */
	public AbnormalCompletionException() {
		this("Completion Abnormal");
	}

	/**
	 * Constructor for MethodNotImplementedException.
	 * @param message
	 */
	public AbnormalCompletionException(String message) {
		super(message);
	}

}
