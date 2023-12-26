package etomica.modules.glass;

import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.IDataSink;
import etomica.data.meter.MeterStructureFactor;

/**
 * This takes structure factor as data and extracts the sin/cos components
 * (using the phase angles) and sends them to a DataSourceCorrelation.
 */
public class StructorFactorComponentExtractor implements IDataSink {
    private double[][] xyData;
    private final MeterStructureFactor meterSFacMobility2;
    private final int interval;
    private final StructureFactorComponentSink dsCorSFacDensityMobility;

    public StructorFactorComponentExtractor(MeterStructureFactor meterSFacMobility2, int interval, StructureFactorComponentSink dsCorSFacDensityMobility) {
        this.meterSFacMobility2 = meterSFacMobility2;
        this.interval = interval;
        this.dsCorSFacDensityMobility = dsCorSFacDensityMobility;
    }

    @Override
    public void putData(IData data) {
        double[] phaseAngles = meterSFacMobility2.getPhaseAngles();
        for (int i = 0; i < xyData.length; i++) {
            double sfac = data.getValue(i);
            double tanphi = Math.tan(phaseAngles[i]);
            double x = Math.sqrt(sfac / (1 + tanphi * tanphi));
            if (phaseAngles[i] > Math.PI / 2 || phaseAngles[i] < -Math.PI / 2) {
                x = -x;
            }
            double y = x * tanphi;
            xyData[i][0] = x;
            xyData[i][1] = y;
        }
        // DataSourceCorrelation expects 1
        dsCorSFacDensityMobility.putData(1, interval, xyData);
    }

    @Override
    public void putDataInfo(IDataInfo dataInfo) {
        xyData = new double[dataInfo.getLength()][2];
    }

    public interface StructureFactorComponentSink {
        // so ugly.  we take idx here because DataSourceCorrelation wants it
        void putData(int idx, int interval, double[][] xyData);
    }
}
