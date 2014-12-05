/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.mvc;

public enum ViewStatus {

  CONTINUE_HELP,
  CONTINUE_PRIOR,
  CONTINUE_NEXT,
  COMPLETE_SUCCESS,
  COMPLETE_USER_CANCELED,
  COMPLETE_EXCEPTION;

  public boolean isEarlyTerminated() {

    return (this == COMPLETE_EXCEPTION) || (this == COMPLETE_USER_CANCELED);
  }

  public boolean isSuccess() {

    return (this == ViewStatus.COMPLETE_SUCCESS);
  }

  public boolean isTerminated() {

    return isSuccess() || isEarlyTerminated();
  }
}