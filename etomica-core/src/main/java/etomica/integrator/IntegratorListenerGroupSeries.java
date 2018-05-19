/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.integrator;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * A composition of IntegratorListeners that must be performed together and in sequence.
 */
public final class IntegratorListenerGroupSeries implements IntegratorListener {

    private final List<IntegratorListener> listeners = new ArrayList<>();
    private int interval;
    private int intervalCount;

    public IntegratorListenerGroupSeries() {
    }

    public IntegratorListenerGroupSeries(IntegratorListener[] listeners) {
        this.listeners.addAll(Arrays.asList(listeners));
        intervalCount = 0;
        interval = 1;
    }

    public void integratorInitialized(IntegratorEvent e) {
        for (IntegratorListener listener : this.listeners) {
            listener.integratorInitialized(e);
        }
    }

    public void integratorStepStarted(IntegratorEvent e) {
        intervalCount++;
        if (intervalCount >= interval) {
            for (IntegratorListener listener : this.listeners) {
                listener.integratorStepStarted(e);
            }
        }
    }

    public void integratorStepFinished(IntegratorEvent e) {
        if (intervalCount >= interval) {
            for (IntegratorListener listener : this.listeners) {
                listener.integratorStepFinished(e);
            }
            intervalCount = 0;
        }
    }

    /**
     * Adds the given listener to the group.  No check is made of whether listener is already in group; it is added
     * regardless.
     *
     * @param newListener
     */
    public void addListener(IntegratorListener newListener) {
        this.listeners.add(newListener);
    }

    /**
     * Removes the given listener from the group.
     */
    public boolean removeListener(IntegratorListener oldListener) {
        return this.listeners.remove(oldListener);
    }

    /**
     * @return an unmodifiable view of the list of listeners.
     */
    public List<IntegratorListener> getListeners() {
        return Collections.unmodifiableList(this.listeners);
    }

    public void setInterval(int i) {
        interval = i;
    }
}
