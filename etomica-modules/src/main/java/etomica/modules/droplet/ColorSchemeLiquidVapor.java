/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.droplet;

import java.awt.Color;

import etomica.atom.IAtom;
import etomica.atom.AtomFilterCollective;
import etomica.graphics.ColorSchemeCollective;

public class ColorSchemeLiquidVapor extends etomica.graphics.ColorScheme implements ColorSchemeCollective {
    
    public ColorSchemeLiquidVapor(AtomFilterCollective liquidFilter) {
        this.liquidFilter = liquidFilter;
        setVaporColor(Color.BLUE);
        setLiquidColor(Color.RED);
    }

    public void setLiquidColor(Color newLiquidColor) {
        liquidColor = newLiquidColor;
    }
    
    public Color getLiquidColor() {
        return liquidColor;
    }

    public void setVaporColor(Color newVaporColor) {
        vaporColor = newVaporColor;
    }
    
    public Color getVaporColor() {
        return vaporColor;
    }
    
    public void setDoResetFilter(boolean newDoResetFilter) {
        doResetFilter = newDoResetFilter;
    }

    public void colorAllAtoms() {
		if (doResetFilter) {
		    liquidFilter.resetFilter();
		}
    }
    
    public Color getAtomColor(IAtom a) {
        return liquidFilter.test(a) ? liquidColor : vaporColor;
    }
    
    protected Color vaporColor, liquidColor;
    protected final AtomFilterCollective liquidFilter;
    protected boolean doResetFilter;
}
