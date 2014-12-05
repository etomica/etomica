/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.mvc;


public enum ActionStatus {

  CONTINUE_SUCCESS,
  COMPLETE_SUCCESS,
  COMPLETE_EXCEPTION;

  public boolean isFailure() {

    return (this == COMPLETE_EXCEPTION);
  }

  public boolean isSuccess() {

    return (this == CONTINUE_SUCCESS) || (this == COMPLETE_SUCCESS);
  }

  public boolean isTerminated() {

    return (this == COMPLETE_EXCEPTION) || (this == COMPLETE_SUCCESS);
  }
}