/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graphics;

import etomica.atom.IAtom;
import etomica.atom.IAtomKinetic;
import etomica.integrator.IntegratorHardFasterer;
import etomica.space.Vector;

import java.awt.*;

/**
 * This colorScheme acts to color differently the two atoms that are scheduled to collide next.
 * Highlight colors are specified by the colliderColor and partnerColor fields; all other
 * atoms are colored with the baseColor.  Applies only to with a hard-potential MD integrator.
 */
public class ColorSchemeColliders extends ColorScheme implements IntegratorHardFasterer.CollisionListener {
    
    IAtom atom1, atom2;

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
    @Override
    public Color getAtomColor(IAtom a) {
        if(a == atom1) return colliderColor;
        if(a == atom2) return partnerColor;
        else return defaultColor;
    }

    @Override
    public void pairCollision(IAtomKinetic atom1, IAtomKinetic atom2, Vector rij, Vector dv, double virial, double tCollision) {
        this.atom1 = atom1;
        this.atom2 = atom2;
    }
}//end of HighlightColliders

