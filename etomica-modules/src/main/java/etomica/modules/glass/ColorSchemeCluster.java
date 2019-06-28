/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.modules.glass;

import etomica.atom.AtomTest;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.graphics.ColorScheme;
import etomica.graphics.ColorSchemeCollective;

import java.awt.*;

public class ColorSchemeCluster extends ColorScheme implements ColorSchemeCollective {
    protected final Box box;
    protected final AtomNbrClusterer clusterer;
    protected final Color[] colors;
    protected final int[][] clusterSize;
    protected final int[] clusterSizeOrder;

    public ColorSchemeCluster(Box box, AtomTest atomTest) {
        this.box = box;
        clusterer = new AtomNbrClusterer(box, atomTest);
        colors = new Color[]{Color.RED, Color.BLUE, Color.GREEN, Color.YELLOW, Color.CYAN, Color.MAGENTA};
        clusterSize = new int[box.getLeafList().size()][2];
        clusterSizeOrder = new int[box.getLeafList().size()];
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
        java.util.Arrays.sort(clusterSize, 0, nClusters, new java.util.Comparator<int[]>() {
            public int compare(int[] a, int[] b) {
                return Integer.compare(a[1], b[1]);
            }
        });
        for (int j = 0; j < nClusters; j++) {
            clusterSizeOrder[clusterSize[nClusters - j - 1][0]] = j;
        }
    }

    @Override
    public Color getAtomColor(IAtom a) {
        int c = clusterSizeOrder[clusterer.getClusters()[a.getLeafIndex()]];
        if (c >= colors.length) return box.getSpace().D() == 2 ? Color.BLACK : Color.WHITE;
        return colors[c];
    }
}
