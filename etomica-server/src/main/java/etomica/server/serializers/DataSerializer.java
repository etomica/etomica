package etomica.server.serializers;

import com.fasterxml.jackson.core.JsonGenerator;
import com.fasterxml.jackson.databind.SerializerProvider;
import com.fasterxml.jackson.databind.ser.std.StdSerializer;
import etomica.data.IData;

import java.io.IOException;

public class DataSerializer extends StdSerializer<IData> {

    public DataSerializer() {
        super(IData.class);
    }

    @Override
    public void serialize(IData value, JsonGenerator gen, SerializerProvider provider) throws IOException {
        double[] dataValues = new double[value.getLength()];
        value.assignTo(dataValues);
        gen.writeArray(dataValues, 0, dataValues.length);
    }
}
