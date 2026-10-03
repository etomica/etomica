/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.modules.glass;

import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.IDataSource;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataTable;
import etomica.integrator.IntegratorMD;
import etomica.space.Vector;
import etomica.units.dimensions.Quantity;
import etomica.util.collections.IntArrayList;
import etomica.util.voro.CLoopAll;
import etomica.util.voro.ContainerPoly;
import etomica.util.voro.PreContainerPoly;
import etomica.util.voro.VoronoiCell;
import etomica.virial.IntSet;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;

public class VoronoiFaceOrders {

    protected final Box box;
    protected final double[] radii;
    protected ContainerPoly con;

    public VoronoiFaceOrders(Box box, double[] radii) {
        this.box = box;
        this.radii = radii;
    }

    private void update() {
        Vector b = box.getBoundary().getBoxSize();
        double bx = 0.5 * b.getX(0);
        double by = 0.5 * b.getX(1);
        double bz = 0.5 * b.getX(2);

        PreContainerPoly pconp = new PreContainerPoly(-bx, bx, -by, by, -bz, bz, true, true, true);
        for (IAtom a : box.getLeafList()) {
            Vector r = a.getPosition();
            pconp.put(a.getLeafIndex(), r.getX(0), r.getX(1), r.getX(2), radii[a.getType().getIndex()]);
        }

        int[] nxout = new int[1];
        int[] nyout = new int[1];
        int[] nzout = new int[1];
        pconp.guess_optimal(nxout, nyout, nzout);
        int nx = nxout[0];
        int ny = nyout[0];
        int nz = nzout[0];

        con = new ContainerPoly(-bx, bx, -by, by, -bz, bz, nx, ny, nz, true, true, true, 8);
        pconp.setup(con);
    }

    public static class DataSourceVoronoiFaceOrders implements IDataSource {

        protected final VoronoiFaceOrders vfo;
        protected final IntegratorMD integrator;
        protected final AtomType atomType;
        protected DataTable data;
        protected DataTable.DataInfoTable dataInfo;
        protected DataTag tag;
        protected final DataDoubleArray.DataInfoDoubleArray[] columnDataInfo;
        protected boolean sortByCount;
        protected long lastStep = -1;

        public DataSourceVoronoiFaceOrders(VoronoiFaceOrders vfo, IntegratorMD integrator, AtomType atomType) {
            this.vfo = vfo;
            this.integrator = integrator;
            this.atomType = atomType;
            data = new DataTable(5, 0);
            columnDataInfo = new DataDoubleArray.DataInfoDoubleArray[5];
            for (int i=0; i<4; i++) {
                columnDataInfo[i] = new DataDoubleArray.DataInfoDoubleArray(""+(i+3), Quantity.DIMENSION, new int[]{0});
            }
            columnDataInfo[4] = new DataDoubleArray.DataInfoDoubleArray("count", Quantity.DIMENSION, new int[]{0});
            dataInfo = new DataTable.DataInfoTable("face orders", columnDataInfo, 0, new String[0]);
            tag = new DataTag();
            dataInfo.addTag(tag);
            sortByCount = true;
        }

         public IData getData() {
            if (integrator.getStepCount() != lastStep) {
                vfo.update();
                lastStep = integrator.getStepCount();
            }
            CLoopAll vl = new CLoopAll(vfo.con);
            HashMap<IntSet,Integer> counts = new HashMap<>();
            if (vl.start()) {
                VoronoiCell c = new VoronoiCell(vfo.con);
                do {
                    int pid = vl.pid();
                    if (atomType != null && atomType != integrator.getBox().getLeafList().get(pid).getType()) continue;
                    if (vfo.con.compute_cell(c, vl)) {
                        IntArrayList v = new IntArrayList();
                        c.face_freq_table(v);
                        int[] ints = new int[4];
                        for (int i=3; i<v.size() && i<=6; i++) {
                            ints[i-3] = v.getInt(i);
                        }
                        IntSet is = new IntSet(ints);
                        counts.merge(is, 1, Integer::sum);
                    }
                } while(vl.inc());
            }
            ArrayList<IntSet> fo = new ArrayList<>(counts.keySet());
            if (sortByCount) {
                fo.sort(new Comparator<IntSet>() {
                    @Override
                    public int compare(IntSet o1, IntSet o2) {
                        int c1 = counts.get(o1);
                        int c2 = counts.get(o2);
                        int rv =  -Integer.compare(c1, c2);
                        if (rv != 0) return rv;
                        return o1.compareTo(o2);
                    }
                });
            }
            else {
                fo.sort(new Comparator<IntSet>() {
                    @Override
                    public int compare(IntSet o1, IntSet o2) {
                        return o1.compareTo(o2);
                    }
                });
            }
            int oldSize = dataInfo.getNRows();
            if (fo.size() > oldSize) {
                data = new DataTable(5, fo.size());
                dataInfo = new DataTable.DataInfoTable("face orders", columnDataInfo, fo.size(), new String[fo.size()]);
                dataInfo.addTag(tag);
            }
            double[][] y = new double[5][];
            for (int i=0; i<5; i++) {
                y[i] = ((DataDoubleArray)data.getData(i)).getData();
            }
            for (int i=0; i<fo.size(); i++) {
                IntSet is = fo.get(i);
                for (int j=0; j<4; j++) y[j][i] = is.v[j];
                y[4][i] = counts.get(is);
            }
            return data;
        }

        @Override
        public DataTag getTag() {
            return tag;
        }

        @Override
        public IDataInfo getDataInfo() {
            return dataInfo;
        }

    }
}
