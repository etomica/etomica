/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.exception;

/**
 * Exception that indicates a method is defined for a class (perhaps to adhere
 * to an interface), but its implementation is not yet completed.
 * 
 * @author David Kofke
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
