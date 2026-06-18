/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.droplet;

import etomica.atom.AtomTestCollective;
import etomica.atom.IAtom;
import etomica.graphics.ColorSchemeCollective;

import java.awt.*;

public class ColorSchemeDropletSurface extends etomica.graphics.ColorScheme implements ColorSchemeCollective {

    public ColorSchemeDropletSurface(AtomTestCollective liquidFilter) {
        this.liquidFilter = liquidFilter;
        setSurfaceColor(Color.RED);
        setInteriorColor(Color.GREEN);
    }

    public void setSurfaceColor(Color newSurfaceColor) {
        surfaceColor = newSurfaceColor;
    }
    
    public Color getSurfaceColor() {
        return surfaceColor;
    }

    public void setInteriorColor(Color newInteriorColor) {
        interiorColor = newInteriorColor;
    }
    
    public Color getInteriorColor() {
        return interiorColor;
    }
    
    public void colorAllAtoms() {
        liquidFilter.resetTest();
    }
    
    public Color getAtomColor(IAtom a) {
        return liquidFilter.test(a) ? interiorColor : surfaceColor;
    }
    
    protected Color surfaceColor, interiorColor;
    protected final AtomTestCollective liquidFilter;
}
