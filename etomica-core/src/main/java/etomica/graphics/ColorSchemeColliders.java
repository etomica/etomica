/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graphics;

import etomica.atom.IAtom;
import etomica.atom.IAtomKinetic;
import etomica.integrator.IntegratorHard;
import etomica.integrator.IntegratorListener;
import etomica.space.Vector;

import java.awt.*;

/**
 * This colorScheme acts to color differently the two atoms that are scheduled to collide next.
 * Highlight colors are specified by the colliderColor and partnerColor fields; all other
 * atoms are colored with the baseColor.  Applies only to with a hard-potential MD integrator.
 */
public class ColorSchemeColliders extends ColorScheme {
    
    IAtom atom1, atom2;

    /**
     * Color applied to the downList atom of the colliding pair
     */
    public java.awt.Color colliderColor = java.awt.Color.red;
    /**
     * Color applied to the upList atom of the colliding pair
     */
    public java.awt.Color partnerColor = java.awt.Color.blue;

    private IntegratorHard integrator;

    public ColorSchemeColliders(IntegratorHard integrator) {
        super();
        this.integrator = integrator;
    }
    /**
     * Applies the special colors to the colliding pair while coloring all other atoms with baseColor.
     */
    @Override
    public Color getAtomColor(IAtom a) {
        IAtom[] colliders = integrator.getNextCollidingPair();
        if(a == colliders[0]) return colliderColor;
        if(a == colliders[1]) return partnerColor;
        else return defaultColor;
    }

}//end of HighlightColliders

