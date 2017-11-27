package etomica.server.representations;

import com.fasterxml.jackson.annotation.JsonProperty;
import etomica.data.IData;
import etomica.data.IDataInfo;

public class DataAndInfo {
    private @JsonProperty
    IData data;
    private @JsonProperty
    IDataInfo dataInfo;

    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    public void setDataInfo(IDataInfo dataInfo) {
        this.dataInfo = dataInfo;
    }

    public IData getData() {
        return data;
    }

    public void setData(IData data) {
        this.data = data;
    }
}
