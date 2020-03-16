/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.modules.glass;

import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.IDataSink;
import etomica.space.Vector;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectOutputStream;

public class DataClusterWriter implements IDataSink {

    protected int nDataP, nSamplesP;
    protected double pSum;
    protected final ObjectOutputStream oos;
    protected final Box box;
    protected float[] x, xyz;

    public DataClusterWriter(Box box, String filename) {
        this.box = box;
        try {
            FileOutputStream fos = new FileOutputStream(filename);
            oos = new ObjectOutputStream(fos);
        } catch (IOException ex) {
            throw new RuntimeException(ex);
        }
    }

    @Override
    public void putData(IData data) {
        if (x.length != data.getLength()) {
            throw new RuntimeException("oops");
        }
        for (int i = 0; i < x.length; i++) {
            x[i] = (float) data.getValue(i);
        }
        float P = (float) (pSum / nSamplesP);
        pSum = nSamplesP = 0;

        if (xyz.length != 3 * box.getLeafList().size()) {
            throw new RuntimeException("oops");
        }
        final Vector p = box.getSpace().makeVector();
        for (IAtom a : box.getLeafList()) {
            int idx = a.getLeafIndex();
            p.E(a.getPosition());
            box.getBoundary().nearestImage(p);
            for (int k = 0; k < 3; k++) {
                xyz[idx * 3 + k] = (float) p.getX(k);
            }
        }

        try {
            oos.writeObject(xyz);
            oos.writeFloat(P);
            oos.writeObject(x);
            // reset allows oos to release the objects (otherwise it will hold on to
            // them so it can write a reference to the previously-written object later)
            oos.reset();
        } catch (IOException ex) {
            throw new RuntimeException(ex);
        }
    }

    public void reset() {
        nDataP = nSamplesP = 0;
        pSum = 0;
    }

    @Override
    public void putDataInfo(IDataInfo inputDataInfo) {
        x = new float[inputDataInfo.getLength()];

        int n = box.getLeafList().size();
        xyz = new float[n * 3];
    }

    public void closeFile() {
        try {
            oos.writeObject(null);
            oos.close();
        } catch (IOException ex) {
            throw new RuntimeException(ex);
        }
    }

    public PressureSink makePressureSink() {
        return new PressureSink();
    }

    public class PressureSink implements IDataSink {
        public void putDataInfo(IDataInfo inputDataInfo) {
        }

        public void putData(IData data) {
            pSum += data.getValue(0);
            nSamplesP++;
        }
    }
}
