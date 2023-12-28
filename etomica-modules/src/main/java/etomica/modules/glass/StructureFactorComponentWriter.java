/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.modules.glass;

import etomica.box.Box;
import etomica.data.ConfigurationStorage;
import etomica.data.meter.MeterStructureFactor;
import etomica.space.Vector;

import java.io.FileWriter;
import java.io.IOException;

public class StructureFactorComponentWriter implements DataSinkBlockAveragerSFac.Sink, StructorFactorComponentExtractor.StructureFactorComponentSink {

    protected final FileWriter fw;
    protected final ConfigurationStorage configStorage;
    protected boolean firstData = true;
    protected final int minInterval;

    public StructureFactorComponentWriter(String filename, MeterStructureFactor meterSFac, ConfigurationStorage configStorage, int minInterval) {
        try {
            fw = new FileWriter(filename);
        }
        catch (IOException ex) {
            throw new RuntimeException(ex);
        }
        this.configStorage = configStorage;
        this.minInterval = minInterval;
        try {
            fw.write("{\"WV\": [");
            Vector L = configStorage.getBox().getBoundary().getBoxSize();
            Vector[] wv = meterSFac.getWaveVectors();
            Vector wvL = configStorage.getBox().getSpace().makeVector();
            for (int j=0; j<wv.length; j++) {
                wvL.E(wv[j]);
                wvL.TE(L);
                wvL.TE(1/(2*Math.PI));
                if (j>0) fw.write(",");
                fw.write(" [");
                boolean firstD = true;
                for (int i=0; i<wvL.getD(); i++) {
                    if (!firstD) fw.write(", ");
                    int iwv = (int)Math.round(wvL.getX(i));
                    fw.write(" "+iwv);
                    firstD = false;
                }
                fw.write(" ]");
            }
            fw.write(" ],\n");
            fw.write("\"data\": [");
        }
        catch (IOException ex) {
            throw new RuntimeException(ex);
        }
    }
    @Override
    public void putData(int interval, double[][] xy) {
        if (interval < minInterval) return;
        long step = configStorage.getSavedSteps()[0];
        try {
            if (!firstData) fw.write(",");
            firstData = false;
            Box box = configStorage.getBox();
            int N = box.getLeafList().size();
            fw.write("{\"step\": " + step + ", \"interval\": " + interval + ", \"sfac\": [");
            for (int i=0; i<xy.length; i++) {
                if (i>0) fw.write(",");
                double sfac = 2*Math.sqrt((xy[i][0]*xy[i][0] + xy[i][1]*xy[i][1])*N);
                double theta = Math.atan2(xy[i][0], xy[i][1]);
                fw.write(String.format("[ %8.3e, %5.3f ]", sfac, theta));
            }
            fw.write("]}\n");
        }
        catch (IOException ex) {
            throw new RuntimeException(ex);
        }

    }

    public void closeFile() {
        try {
            fw.write("]}");
            fw.close();
        }
        catch (IOException ex) {
            throw new RuntimeException(ex);
        }
    }

    @Override
    public void putData(int idx, int interval, double[][] xyData) {
        // so ugly.  we take idx here because DataSourceCorrelation wants it
        putData(interval, xyData);
    }
}
