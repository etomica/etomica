/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.modules.nucleation;

import etomica.atom.AtomTest;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.graphics.ColorScheme;
import etomica.graphics.ColorSchemeCollective;
import etomica.modules.glass.AtomNbrClusterer;

import java.awt.*;

public class ColorSchemeCluster extends ColorScheme implements ColorSchemeCollective {
    protected final Box box;
    protected final AtomNbrClusterer clusterer;
    protected final int[][] clusterSize;

    public ColorSchemeCluster(Box box, AtomTest atomTest) {
        this.box = box;
        clusterer = new AtomNbrClusterer(box, atomTest);
        clusterSize = new int[box.getLeafList().size()][2];
    }

    public void setNbrMax(double nbrMax) {
        clusterer.setNbrMax(nbrMax);
    }

    public double getNbrMax() {
        return clusterer.getNbrMax();
    }

    @Override
    public void colorAllAtoms() {
        clusterer.findClusters();
        int[] firstAtom = clusterer.getFirstAtom();
        int[] nextAtom = clusterer.getNextAtom();
        int nClusters = 0;
        for (int i = 0; i < firstAtom.length; i++) {
            if (firstAtom[i] == -1) break;
            nClusters++;
            clusterSize[i][0] = i;
            clusterSize[i][1] = 0;
            for (int ii = firstAtom[i]; ii != -1; ii = nextAtom[ii]) {
                clusterSize[i][1]++;
            }
        }
    }

    @Override
    public Color getAtomColor(IAtom a) {
        int c = clusterer.getClusters()[a.getLeafIndex()];
        int s = clusterSize[c][1];
        int n = box.getLeafList().size();
        float x = (float) Math.sqrt(Math.log(s) / Math.log(n));
        if (x < 0.5) return new Color(1, 0, 2 * x);
        else return new Color(2 * (1 - x), 0, 1);
    }
}
