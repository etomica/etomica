package etomica.modules.droplet;

import java.awt.Color;

import etomica.api.IAtomLeaf;
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
    
    public Color getAtomColor(IAtomLeaf a) {
        return liquidFilter.accept(a) ? liquidColor : vaporColor;
    }
    
    private static final long serialVersionUID = 1L;
    protected Color vaporColor, liquidColor;
    protected final AtomFilterCollective liquidFilter;
    protected boolean doResetFilter;
}
 