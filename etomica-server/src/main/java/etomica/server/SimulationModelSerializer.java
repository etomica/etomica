package etomica.server;

import com.fasterxml.jackson.core.JsonGenerator;
import com.fasterxml.jackson.databind.SerializerProvider;
import com.fasterxml.jackson.databind.ser.std.StdSerializer;
import etomica.meta.SimulationModel;
import etomica.meta.wrappers.Wrapper;

import java.io.IOException;
import java.util.function.Function;
import java.util.stream.Collectors;

public class SimulationModelSerializer extends StdSerializer<SimulationModel> {


    protected SimulationModelSerializer(Class<SimulationModel> t) {
        super(t);
    }

    @Override
    public void serialize(SimulationModel value, JsonGenerator gen, SerializerProvider provider) throws IOException {
        gen.writeObject(
                value.allWrappers().stream().collect(Collectors.toMap(Wrapper::getWrappedId, Function.identity()))
        );

    }
}
