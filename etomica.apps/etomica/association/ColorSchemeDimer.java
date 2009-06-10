package etomica.association;

import java.awt.Color;

import etomica.api.IAtom;
import etomica.api.IBox;
import etomica.api.IRandom;
import etomica.atom.AtomLeafAgentManager;
import etomica.graphics.ColorScheme;

/**
 * Color scheme for dimer association.  All monomers have the same color and
 * Each dimer has a different color.
 * 
 * @author Andrew Schultz
 */
public class ColorSchemeDimer extends ColorScheme implements AtomLeafAgentManager.AgentSource {

    public ColorSchemeDimer(AssociationManager associationManager, IBox box, IRandom random) {
        super();
        setMonomerColor(DEFAULT_ATOM_COLOR);
        this.associationManager = associationManager;
        dimerColorManager = new AtomLeafAgentManager(this, box);
        this.random = random;
    }
    
    public Color getAtomColor(IAtom a) {
        if (associationManager.getAssociatedAtoms(a).getAtomCount() == 0) {
            return monomerColor;
        }
        Color color1 = (Color)dimerColorManager.getAgent(a);
        if (color1 != null && color1 == dimerColorManager.getAgent(associationManager.getAssociatedAtoms(a).getAtom(0))) {
            return color1;
        }
        Color newColor = new Color((float)random.nextDouble(),(float)random.nextDouble(),(float)random.nextDouble());
        dimerColorManager.setAgent(a, newColor);
        dimerColorManager.setAgent(associationManager.getAssociatedAtoms(a).getAtom(0), newColor);
        return newColor;
    }
    
    /**
     * Sets the color used for all dimers.
     */
    public void setMonomerColor(Color newMonomerColor) {
        monomerColor = newMonomerColor;
    }

    /**
     * Returns the color used for all dimers.
     */
    public Color getMonomerColor() {
        return monomerColor;
    }
    
    public Class getAgentClass() {
        return Color.class;
    }

    public Object makeAgent(IAtom a) {
        return null;
    }

    public void releaseAgent(Object agent, IAtom atom) {
    }
    
    private static final long serialVersionUID = 1L;
    protected Color dimerColor, monomerColor;
    protected AssociationManager associationManager;
    protected AtomLeafAgentManager dimerColorManager;
    protected IRandom random;
}
