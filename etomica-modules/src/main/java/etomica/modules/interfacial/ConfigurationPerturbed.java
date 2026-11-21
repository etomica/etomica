/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.modules.interfacial;

import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.config.Configuration;
import etomica.space.Vector;
import etomica.util.random.IRandom;

public class ConfigurationPerturbed implements Configuration {

    protected final Configuration config;
    protected final IRandom random;
    protected double d;

    public ConfigurationPerturbed(Configuration config, IRandom random, double d) {
        this.config = config;
        this.random = random;
        this.d = d;
    }

    @Override
    public void initializeCoordinates(Box box) {
        config.initializeCoordinates(box);
        Vector dr = box.getSpace().makeVector();
        for (IAtom a : box.getLeafList()) {
            dr.setRandomInSphere(random);
            a.getPosition().PEa1Tv1(d, dr);
        }
    }
}
