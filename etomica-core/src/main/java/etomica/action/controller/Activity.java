package etomica.action.controller;

public abstract class Activity {

    public void restart() { }

    public abstract void runActivity(Controller.ControllerHandle handle);

}
