package etomica.server;

import com.fasterxml.jackson.core.JsonGenerator;
import com.fasterxml.jackson.databind.SerializerProvider;
import com.fasterxml.jackson.databind.ser.std.StdSerializer;
import etomica.meta.SimulationModel;

import java.io.IOException;

public class SimulationModelSerializer extends StdSerializer<SimulationModel> {


    protected SimulationModelSerializer(Class<SimulationModel> t) {
        super(t);
    }

    @Override
    public void serialize(SimulationModel value, JsonGenerator gen, SerializerProvider provider) throws IOException {
        gen.writeObject(value.allWrappers());

    }
}
