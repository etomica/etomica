/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.association;

import java.awt.Color;
import java.util.HashMap;

import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.util.random.IRandom;
import etomica.association.ColorSchemeSmer.ColorAgent;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomLeafAgentManager;
import etomica.graphics.ColorScheme;

/**
 * Color scheme for smer association.  All monomers have the same color and
 * Each smer has a different color.
 * 
 * @author Andrew Schultz
 */
public class ColorSchemeSmer extends ColorScheme implements AtomLeafAgentManager.AgentSource<ColorAgent> {

    public ColorSchemeSmer(IAssociationHelper associationHelper, Box box, IRandom random) {
        super();
        setMonomerColor(DEFAULT_ATOM_COLOR);
        this.associationHelper = associationHelper;
        dimerColorManager = new AtomLeafAgentManager<ColorAgent>(this, box, ColorAgent.class);
        this.random = random;
        oldColors = new HashMap<Color,Integer>();
        smerList = new AtomArrayList();
    }
    
    protected boolean equalLists(AtomArrayList list1, AtomArrayList list2) {
        for (int i=0; i<list1.getAtomCount(); i++) {
            if (list1.indexOf(list2.getAtom(i)) == -1) {
                return false;
            }
        }
        return true;
    }

    protected void reclaim(IAtom atom) {
//        System.out.println("reclaiming "+atom);
        ColorAgent colorAgent = dimerColorManager.getAgent(atom);
        if (colorAgent == null) return;
        dimerColorManager.setAgent(atom, null);
        if (oldColors.get(colorAgent.color) == null) {
//            System.out.println("returning "+colorAgent.color+" to the hash");
            oldColors.put(colorAgent.color, 1);
        }
        for (int i=0; i<colorAgent.smerList.getAtomCount(); i++) {
            if (colorAgent.smerList.getAtom(i) == atom) continue;
            reclaim(colorAgent.smerList.getAtom(i));
        }
    }
    
    public synchronized Color getAtomColor(IAtom a) {
        ColorAgent colorAgent = dimerColorManager.getAgent(a);
        associationHelper.populateList(smerList, a, false);
        if (colorAgent != null) {
            if (smerList.getAtomCount() > 1) {
                if (smerList.getAtomCount() != colorAgent.smerList.getAtomCount() || !equalLists(smerList, colorAgent.smerList)) {
                    // atom was in an smer, and is now in a different smer (perhaps shrunk, grown, combined, etc)
                    // discard info for old smer
//                    System.out.println("smer gone "+colorAgent.smerList);
                    // reclaim a and all atoms that used to be with a
                    reclaim(a);
                    
//                    System.out.println(" => new smer "+smerList+"  reusing "+colorAgent.color);
                    // use colorAgent for new smer.  keep the old color
                    for (int i=1; i<smerList.getAtomCount(); i++) {
                        reclaim(smerList.getAtom(i));
                    }
                    colorAgent.smerList.clear();
                    colorAgent.smerList.addAll(smerList);
                    // the act of reclaiming returned our color to the hash.  retrieve a new one
                    colorAgent.color = (Color)oldColors.keySet().toArray()[0];
                    oldColors.remove(colorAgent.color);
                    for (int i=1; i<smerList.getAtomCount(); i++) {
                        dimerColorManager.setAgent(colorAgent.smerList.getAtom(i), colorAgent);
                    }
                        
                    return colorAgent.color;
                }
                // atom is in the same smer as before
                return colorAgent.color;
            }
            // atom was in an smer, and is now a monomer
//            System.out.println("smer totally gone "+colorAgent.smerList);
            // discard color and info from old smerList
            reclaim(a);
            
            return monomerColor;
        }
        // atom was a monomer
        if (smerList.getAtomCount() == 1) {
            // atom is still a monomer
            return monomerColor;
        }
        // atom is in an smer.  see if any of the atoms were in an smer
        ColorAgent iColorAgent = null;
        for (int i=1; i<smerList.getAtomCount(); i++) {
            iColorAgent = dimerColorManager.getAgent(smerList.getAtom(i));
            if (iColorAgent != null) {
                colorAgent = iColorAgent;
                break;
            }
        }
        if (colorAgent != null) {
            // discard info for old smer
//            System.out.println("smer gone "+colorAgent.smerList);
            for (int i=0; i<colorAgent.smerList.getAtomCount(); i++) {
                reclaim(colorAgent.smerList.getAtom(i));
            }
            
            // use colorAgent for new smer.  keep the old color
//            System.out.println(" => new smer "+smerList+"  reusing "+colorAgent.color);
            colorAgent.smerList.clear();
            colorAgent.smerList.addAll(smerList);
            // the act of reclaiming returned our color to the hash.  retrieve a new one
            colorAgent.color = (Color)oldColors.keySet().toArray()[0];
            oldColors.remove(colorAgent.color);
            for (int i=0; i<smerList.getAtomCount(); i++) {
                dimerColorManager.setAgent(colorAgent.smerList.getAtom(i), colorAgent);
            }

            return colorAgent.color;
        }
//        System.out.println("new new smer "+smerList);
        // none of the new atoms were in an smer
        colorAgent = new ColorAgent();
        // use a random color
        if (oldColors.size() > 0) {
            // retrieve a color from the hash
            colorAgent.color = (Color)oldColors.keySet().toArray()[0];
//            System.out.println("reusing "+colorAgent.color+" from hash");
            oldColors.remove(colorAgent.color);
        }
        else {
            // hash is empty, make a new color
            colorAgent.color = new Color((float)random.nextDouble(), (float)random.nextDouble(), (float)random.nextDouble());

//            System.out.println("made a new color "+colorAgent.color);
        }
        colorAgent.smerList = new AtomArrayList();
        colorAgent.smerList.addAll(smerList);

        // use colorAgent for new smer.
        for (int i=0; i<smerList.getAtomCount(); i++) {
            dimerColorManager.setAgent(colorAgent.smerList.getAtom(i), colorAgent);
        }

        return colorAgent.color;
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

    public ColorAgent makeAgent(IAtom a, Box agentBox) {
        return null;
    }

    public void releaseAgent(ColorAgent agent, IAtom atom, Box agentBox) {
    }
    
    private static final long serialVersionUID = 1L;
    protected Color monomerColor;
    protected IAssociationHelper associationHelper;
    protected AtomLeafAgentManager<ColorAgent> dimerColorManager;
    protected IRandom random;
    protected HashMap<Color,Integer> oldColors;
    protected final AtomArrayList smerList;
    
    public static class ColorAgent {
        public Color color;
        public AtomArrayList smerList;
    }
}
