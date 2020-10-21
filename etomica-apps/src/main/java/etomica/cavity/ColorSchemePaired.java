/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.cavity;

import etomica.atom.IAtom;
import etomica.graphics.ColorScheme;

import java.awt.*;

public class ColorSchemePaired extends ColorScheme {

    protected final P2HardSphereCavity p2;

    public ColorSchemePaired(P2HardSphereCavity p2) {
        this.p2 = p2;
    }

    @Override
    public Color getAtomColor(IAtom a) {
        boolean paired = p2.getPairedAtom1() == a || p2.getPairedAtom2() == a;
        return paired ? Color.RED : Color.BLUE;
    }
}
