/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graphics;
import java.awt.Color;

import etomica.atom.IAtom;
import etomica.integrator.IntegratorHard;

/**
 * This colorScheme acts to color differently the two atoms that are scheduled to collide next.
 * Highlight colors are specified by the colliderColor and partnerColor fields; all other
 * atoms are colored with the baseColor.  Applies only to with a hard-potential MD integrator.
 */
public class ColorSchemeColliders extends ColorScheme {
    
    IntegratorHard integrator;
    
    public ColorSchemeColliders(IntegratorHard integrator) {
        super();
        this.integrator = integrator;
    }
    /**
     * Color applied to the downList atom of the colliding pair
     */
    public java.awt.Color colliderColor = java.awt.Color.red;
    /**
     * Color applied to the upList atom of the colliding pair
     */
    public java.awt.Color partnerColor = java.awt.Color.blue;
    /**
     * Applies the special colors to the colliding pair while coloring all other atoms with baseColor.
     */ 
    public Color getAtomColor(IAtom a) {
        IntegratorHard.Agent colliderAgent = integrator.getLastColliderAgent();
        if(colliderAgent == null) return defaultColor;
        else if(a == colliderAgent.atom) return colliderColor;
        else if(a == colliderAgent.collisionPartner) return partnerColor;
        else return defaultColor;
    }
}//end of HighlightColliders

