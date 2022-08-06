package etomica.modules.glass;

import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.IDataSink;
import etomica.data.meter.MeterStructureFactor;
import etomica.util.Statefull;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class DataSinkBlockAveragerSFac implements IDataSink, Statefull {

    protected final ConfigurationStorage configStorage;
    protected final MeterStructureFactor meterSFac;
    protected final List<Sink> sinks;
    protected double[][][] blockAvg;
    protected int minInterval;
    protected long[] nSamples;
    protected int nData;

    public DataSinkBlockAveragerSFac(ConfigurationStorage configStorage, int minInterval, MeterStructureFactor meterSFac) {
        this.configStorage = configStorage;
        this.meterSFac = meterSFac;
        this.sinks = new ArrayList<>();
        blockAvg = new double[0][0][0];
        this.minInterval = minInterval;
        nSamples = new long[0];
        reallocate(minInterval + 1);
    }

    public void addSink(Sink s) {
        sinks.add(s);
    }

    protected void reallocate(int n) {
        int nOld = blockAvg.length;
        blockAvg = Arrays.copyOf(blockAvg, n);
        for (int i = nOld; i < n; i++) {
            blockAvg[i] = new double[nData][2];
        }
        nSamples = Arrays.copyOf(nSamples, n);
    }

    public void putData(IData inputData) {
        long step = configStorage.getSavedSteps()[0];
        double[][] xy = blockAvg[minInterval];
        double[] phaseAngles = meterSFac.getPhaseAngles();
        for (int j = 0; j < xy.length; j++) {
            double sfac = inputData.getValue(j);
            double tanphi = Math.tan(phaseAngles[j]);
            double x = Math.sqrt(sfac / (1 + tanphi * tanphi));
            if (phaseAngles[j] > Math.PI / 2 || phaseAngles[j] < -Math.PI / 2) {
                x = -x;
            }
            double y = x * tanphi;
            blockAvg[minInterval][j][0] = x;
            blockAvg[minInterval][j][1] = y;
        }
        for (Sink s : sinks) {
            s.putData(minInterval, blockAvg[minInterval]);
        }
        for (int i = minInterval + 1; i <= configStorage.getLastConfigIndex(); i++) {
            if (i >= blockAvg.length) reallocate(i + 1);
            xy = blockAvg[i];
            nSamples[i]++;
            for (int j = 0; j < xy.length; j++) {
                xy[j][0] += (blockAvg[i - 1][j][0] - xy[j][0]) / nSamples[i];
                xy[j][1] += (blockAvg[i - 1][j][1] - xy[j][1]) / nSamples[i];
            }
            if (nSamples[i] == 1) break;
        }
        for (int i = minInterval + 1; i < configStorage.getLastConfigIndex(); i++) {
            if (step % (1L << i) == 0) {
                if (nSamples[i] == 2) {
                    for (Sink s : sinks) {
                        s.putData(i, blockAvg[i]);
                    }
                }
                for (int j = 0; j < blockAvg[i].length; j++) {
                    Arrays.fill(blockAvg[i][j], 0);
                }
                nSamples[i] = 0;
            }
        }
    }

    @Override
    public void putDataInfo(IDataInfo dataInfo) {
        nData = dataInfo.getLength();
        for (int i = 0; i < blockAvg.length; i++) {
            blockAvg[i] = new double[nData][2];
            nSamples[i] = 0;
        }
    }

    @Override
    public void saveState(Writer fw) throws IOException {
        fw.write(nSamples.length+"\n");
        for (int i=0; i<nSamples.length; i++) {
            fw.write(""+nSamples[i]);
            for (int j=0; j<blockAvg[i].length; j++) {
                fw.write(" "+blockAvg[i][j][0]+" "+blockAvg[i][j][1]);
            }
            fw.write("\n");
        }
    }

    @Override
    public void restoreState(BufferedReader br) throws IOException {
        int n = Integer.parseInt(br.readLine());
        reallocate(n);
        for (int i=0; i<n; i++) {
            String[] bits = br.readLine().split(" ");
            nSamples[i] = Long.parseLong(bits[0]);
            for (int j=0; j<blockAvg[i].length; j++) {
                blockAvg[i][j][0] = Double.parseDouble(bits[1+2*j]);
                blockAvg[i][j][1] = Double.parseDouble(bits[1+2*j+1]);
            }
        }
    }

    public interface Sink {
        void putData(int interval, double[][] xy);
    }
}
