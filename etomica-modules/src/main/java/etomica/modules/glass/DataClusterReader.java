package etomica.modules.glass;

import etomica.data.IDataSink;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDoubleArray;
import etomica.units.dimensions.Null;
import etomica.units.dimensions.Pressure;

import java.io.FileInputStream;
import java.io.IOException;
import java.io.ObjectInputStream;

public class DataClusterReader {

    public static void readClusterFile(DataClusterer clusterer, String filename, int skip) {
        try {
            DataDouble.DataInfoDouble pDataInfo = new DataDouble.DataInfoDouble("pressure", Pressure.DIMENSION);
            IDataSink pSink = clusterer.makePressureSink();
            pSink.putDataInfo(pDataInfo);
            DataDouble pData = new DataDouble();

            DataDoubleArray.DataInfoDoubleArray sDataInfo = null;
            DataDoubleArray sData = null;
            double[] s = null;

            FileInputStream fis = new FileInputStream(filename);
            ObjectInputStream in = new ObjectInputStream(fis);

            int countdown = skip;
            while (true) {
                float[] xyz = (float[]) in.readObject();
                if (xyz == null) break;
                pData.x = in.readFloat();
                if (countdown == skip) pSink.putData(pData);
                float[] x = (float[]) in.readObject();
                if (sDataInfo == null) {
                    sDataInfo = new DataDoubleArray.DataInfoDoubleArray("sfac", Null.DIMENSION, new int[]{x.length});
                    clusterer.putDataInfo(sDataInfo);
                    sData = new DataDoubleArray(x.length);
                    s = sData.getData();
                }
                for (int i = 0; i < x.length; i++) {
                    s[i] = x[i];
                }
                if (countdown == skip) clusterer.putData(sData);
                countdown--;
                if (countdown == 0) countdown = skip;
            }
            in.close();
            fis.close();
        } catch (IOException | ClassNotFoundException e) {
            throw new RuntimeException(e);
        }
    }
}
