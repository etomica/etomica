/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graphics;
import java.awt.Color;

import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.box.storage.ObjectStorage;
import etomica.box.storage.Token;
import etomica.box.storage.Tokens;
import etomica.space.Space;
import etomica.util.random.IRandom;

public class ColorSchemeRandom extends ColorScheme {
    
    public ColorSchemeRandom(Box box, IRandom random) {
    	super();
        Token<ObjectStorage<Color>> colors = Tokens.objects(
                new ObjectStorage.Factory<Color>() {
                    @Override
                    public Class<? extends Color> getVClass() {
                        return Color.class;
                    }

                    @Override
                    public Color create(int idx) {
                        return new Color((float) random.nextDouble(), (float) random.nextDouble(), (float) random.nextDouble());
                    }
                }
        );
        colorStorage = box.getAtomStorage(colors);
    }
    
    public Color getAtomColor(IAtom a) {
        return colorStorage.get(a.getLeafIndex());
    }
   
    private final ObjectStorage<Color> colorStorage;
}
