package etomica.association;

import java.awt.Color;
import java.util.HashMap;

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
        oldColors = new HashMap<Color,Integer>();
    }
    
    public Color getAtomColor(IAtom a) {
        Color color1 = (Color)dimerColorManager.getAgent(a);
        if (associationManager.getAssociatedAtoms(a).getAtomCount() == 0) {
            if (color1 != null) {
                if (oldColors.get(color1) == null) {
                    // dimer bond was broken, return color to the hash for re-use
                    oldColors.put(color1, 1);
                }
                // remember that atom is a monomer
                dimerColorManager.setAgent(a, null);
            }
            return monomerColor;
        }
        if (color1 != null && color1 == dimerColorManager.getAgent(associationManager.getAssociatedAtoms(a).getAtom(0))) {
            // we already assigned a color
            return color1;
        }
        if (color1 != null && oldColors.get(color1) == null) {
            // atom formed a new dimer with a different partner, put old color in the hash
            oldColors.put(color1, 1);
        }
        Color color2 = (Color)dimerColorManager.getAgent(associationManager.getAssociatedAtoms(a).getAtom(0));
        if (color2 != null && oldColors.get(color2) == null) {
            // bonding parter had an old color from an old dimer, put old color in the hash
            oldColors.put(color2, 1);
        }
        Color newColor;
        if (oldColors.size() > 0) {
            // retrieve a color from the hash
            newColor = (Color)oldColors.keySet().toArray()[0];
            oldColors.remove(newColor);
        }
        else {
            // hash is empty, make a new color
            newColor = new Color((float)random.nextDouble(),(float)random.nextDouble(),(float)random.nextDouble());
        }
        // remember our new color for bonded atoms
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
    protected HashMap<Color,Integer> oldColors;
}
