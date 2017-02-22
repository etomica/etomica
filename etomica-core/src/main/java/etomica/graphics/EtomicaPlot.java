/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graphics;

import java.awt.Graphics;

import ptolemy.plot.Plot;

/**
 * Custom Plot subclass that allows us to sanely reset the axes autoscaling and
 * repaint when a hidden plot is shown (without resetting all plot parameters).
 * 
 * @author Andrew Schultz
 */
public class EtomicaPlot extends Plot {
    
    public EtomicaPlot(DisplayPlot displayPlot) {
        this.displayPlot = displayPlot;
    }
    
    /**
     * Overrides a not-so-useful method in superclass
     */
    public void resetAxes() {
        // We need to repaint the offscreen buffer.
        _plotImage = null;

        _xBottom = Double.MAX_VALUE;
        _xTop = - Double.MAX_VALUE;
        _yBottom = Double.MAX_VALUE;
        _yTop = - Double.MAX_VALUE;
        if (true) {
            // Protected members first.
            _yMax = 0;
            _yMin = 0;
            _xMax = 0;
            _xMin = 0;
            _xRangeGiven = false;
            _yRangeGiven = false;
            _rangesGivenByZooming = false;
        }
    }

    /**
     * Notify that the displayPlot has new data, but we're not showing.
     * Remember and ask for the data when we are showing.
     */
    public void notifyNeedUpdate() {
        needUpdate = true;
    }

    public synchronized void paintComponent(Graphics graphics) {
        if (needUpdate) {
            // DisplayPlot has updated data for us; the data arrived
            // when we were hidden and now we're shown
            displayPlot.doUpdate();
            needUpdate = false;
        }
        else {
            super.paintComponent(graphics);
        }
    }
    
    public Graphics getGraphics() {
        // this is so dumb, this is the only way to add points to the plot
        // while the plot isn't showing.
        if (!isShowing()) {
            return null;
        }
        return super.getGraphics();
    }
    
    protected final DisplayPlot displayPlot;
    protected boolean needUpdate = false;
}
