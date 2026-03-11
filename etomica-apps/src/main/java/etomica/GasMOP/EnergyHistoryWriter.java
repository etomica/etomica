package etomica.GasMOP;

import etomica.data.*;

import java.io.FileWriter;
import java.io.PrintWriter;

class EnergyHistoryWriter implements IDataSink {


    private final PrintWriter out;
    private final DataSourceCountTime timer;

    public EnergyHistoryWriter(FileWriter writer, DataSourceCountTime timer) {
        this.out = new PrintWriter(writer);
        this.timer = timer;
    }

    // Called once at beginning with metadata
    public void putDataInfo(IDataInfo dataInfo) {
        // Optional: write header
        out.println("# time energy");
    }

    // Called every time energy is pushed
    @Override
    public void putData(IData data) {

        double t = timer.getData().getValue(0);
        double U = data.getValue(0);

        out.println(t + " " + U);
    }

    public DataTag getTag() {
        return null;
    }

    public void close() {
        out.flush();
        out.close();
    }
}
