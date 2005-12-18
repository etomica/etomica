package etomica.graphics;
import java.awt.Color;

import etomica.atom.Atom;
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
    public Color getAtomColor(Atom a) {
        IntegratorHard.Agent colliderAgent = integrator.colliderAgent();
        if(colliderAgent == null) return baseColor;
        else if(a == colliderAgent.atom) return colliderColor;
        else if(a == colliderAgent.collisionPartner) return partnerColor;
        else return baseColor;
    }
}//end of HighlightColliders

