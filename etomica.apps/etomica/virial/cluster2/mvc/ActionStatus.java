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