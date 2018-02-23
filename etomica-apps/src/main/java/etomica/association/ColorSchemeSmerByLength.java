/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.association;

import java.awt.Color;
import java.util.ArrayList;

import etomica.atom.IAtom;
import etomica.util.random.IRandom;
import etomica.atom.AtomArrayList;
import etomica.graphics.ColorScheme;

/**
 * Color scheme for smer association.  All smers of a length have the same
 * color.
 * 
 * @author Andrew Schultz
 */
public class ColorSchemeSmerByLength extends ColorScheme {

    public ColorSchemeSmerByLength(IAssociationHelper associationHelper, IRandom random) {
        super();
        this.associationHelper = associationHelper;
        this.random = random;
        smerList = new AtomArrayList();
        colors = new ArrayList<Color>(2);
        colors.add(null);
        colors.add(DEFAULT_ATOM_COLOR);
        colors.add(new Color(128, 200, 0));
        colors.add(new Color(0, 180, 180));
        colors.add(new Color(128, 0, 200));
        colors.add(new Color(0, 128, 0));
    }
    
    public synchronized Color getAtomColor(IAtom a) {
        associationHelper.populateList(smerList, a, false);
        int oldSize = colors.size();
        if (oldSize < smerList.size()+1) {
            for (int i = oldSize; i<smerList.size()+1; i++) {
                colors.add(new Color((float)random.nextDouble(), (float)random.nextDouble(), (float)random.nextDouble()));
            }
        }
        return colors.get(smerList.size());
    }
    
    private static final long serialVersionUID = 1L;
    protected IAssociationHelper associationHelper;
    protected IRandom random;
    protected final ArrayList<Color> colors;
    protected final AtomArrayList smerList;
}
