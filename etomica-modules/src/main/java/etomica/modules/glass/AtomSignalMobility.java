package etomica.modules.glass;

import etomica.atom.IAtom;
import etomica.data.ConfigurationStorage;
import etomica.data.meter.MeterStructureFactor;
import etomica.space.Vector;

public class AtomSignalMobility extends MeterStructureFactor.AtomSignalSourceByType implements DataSinkBlockAveragerSFac.Sink {
    protected final ConfigurationStorage configStorage;
    protected final Vector dr;
    protected int prevConfigIndex;
    protected final MeterStructureFactor meterDensity;
    protected final int N;
    protected double[][] savedXY;
    protected long savedStep;

    public AtomSignalMobility(ConfigurationStorage configStorage) {
        this(configStorage, null, 0);
    }

    public AtomSignalMobility(ConfigurationStorage configStorage, MeterStructureFactor meterDensity, int N) {
        this.configStorage = configStorage;
        this.meterDensity = meterDensity;
        dr = configStorage.getBox().getSpace().makeVector();
        this.N = N;
        savedXY = new double[0][0];
    }

    public boolean ready() {
        int lastIndex = configStorage.getLastConfigIndex();
        return prevConfigIndex <= lastIndex;
    }

    public void setPrevConfig(int prevConfigIndex) {
        this.prevConfigIndex = prevConfigIndex;
    }

    public int getPrevConfigIndex() {
        return prevConfigIndex;
    }

    public double signal(IAtom atom, int iq) {
        double s = super.signal(atom);
        if (s == 0) return 0;
        int idx = prevConfigIndex;
        int lastIndex = configStorage.getLastConfigIndex();
        if (lastIndex < idx) throw new RuntimeException("not enough configs to compute signal");
        Vector[] positions = configStorage.getSavedConfig(0);
        Vector[] prevPositions = configStorage.getSavedConfig(idx);
        int atomIndex = atom.getLeafIndex();
        dr.Ev1Mv2(positions[atomIndex], prevPositions[atomIndex]);
        double r2 = dr.squared();
        s *= r2;

        if (meterDensity==null) return s;
        if (savedStep != configStorage.getSavedSteps()[0]) throw new RuntimeException("oops "+prevConfigIndex+" wrong step "+savedStep+" "+configStorage.getSavedSteps()[0]);
        double beta = 2*Math.sqrt((savedXY[iq][0]*savedXY[iq][0] + savedXY[iq][1]*savedXY[iq][1])*configStorage.getBox().getLeafList().size());
        double theta = Math.atan2(savedXY[iq][0], savedXY[iq][1]);
        dr.Ev1Pv2(prevPositions[atomIndex], positions[atomIndex]);
        double V = configStorage.getBox().getBoundary().volume();
        Vector wv = meterDensity.getWaveVectors()[iq];
        double den = N/V + beta/V*Math.sin(wv.dot(dr) + theta);
        return s/den;
    }

    @Override
    public void putData(int interval, double[][] xy) {
        if (interval != prevConfigIndex-1) return;
        if (savedXY.length != xy.length) {
            savedXY = new double[xy.length][2];
        }
        for (int i=0; i<xy.length; i++) {
            savedXY[i][0] = xy[i][0];
            savedXY[i][1] = xy[i][1];
        }
        savedStep = configStorage.getSavedSteps()[0];
    }
}
