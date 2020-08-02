package etomica.action.activity;

import etomica.action.IAction;

public interface Activity2 extends IAction {

    void preAction();

    void postAction();

    void restart();
}
