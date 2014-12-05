/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.mvc;

public class MVCException extends Exception {

  private static final long serialVersionUID = -2298180432380683045L;

  public MVCException() {

    super();
  }

  public MVCException(Throwable cause) {

    super(cause);
  }

  public MVCException(String message) {

    super(message);
  }

  public MVCException(String message, Throwable cause) {

    super(message);
  }
}
