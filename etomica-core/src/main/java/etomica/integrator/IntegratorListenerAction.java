/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.integrator;

import etomica.action.IAction;

/**
 * An IntegratorListener that performs an action after every integrator step,
 * optionally with a given interval of integrator steps.
 */
public final class IntegratorListenerAction implements IntegratorListener {

    private final IAction action;
    private int interval;
    private int intervalCount;

    public IntegratorListenerAction(IAction a) {
        this(a, 1);
    }

    public IntegratorListenerAction(IAction a, int interval) {
        action = a;
        this.interval = interval;
        intervalCount = 0;
    }

    public void integratorInitialized(IntegratorEvent e) {
    }

    public void integratorStepStarted(IntegratorEvent e) {
        ++intervalCount;
    }

    public void integratorStepFinished(IntegratorEvent e) {
        if (intervalCount >= interval) {
            intervalCount = 0;
            action.actionPerformed();
        }
    }

    public int getInterval() {
        return interval;
    }

    public void setInterval(int i) {
        interval = i;
    }

    public IAction getAction() {
        return action;
    }

}
