/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.statistics;


import etomica.data.history.HistoryScrolling;

/**
 * History that records a number of values, with new ones replacing the
 * earliest ones when the record is full.  The data returned by the
 * getHistory method will have the earliest first in the array, and the
 * most recent last.  If presented as a plot, the effect will be to
 * scroll the data across the plot window.
 *
 * @author kofke, schultz, cribbin
 */
public class HistoryStatistics extends HistoryScrolling {

    public HistoryStatistics() {
        this(100);
    }

    public HistoryStatistics(int n) {
        super(n);
    }

    public double[] getXValues() {
        int n = history.length;

        for (int i = 0; i < n; i++) {
            tempX[i] = i + 1 - n;
        }

        return tempX;
    }
}
