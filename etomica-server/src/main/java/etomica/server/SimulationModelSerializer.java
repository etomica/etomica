package etomica.server;

import com.fasterxml.jackson.core.JsonGenerator;
import com.fasterxml.jackson.databind.SerializerProvider;
import com.fasterxml.jackson.databind.ser.std.StdSerializer;
import etomica.meta.InstanceProperty;

import java.io.IOException;
import java.util.List;
import java.util.Map;

public class SimulationModelSerializer extends StdSerializer<SimulationModel> {


    protected SimulationModelSerializer(Class<SimulationModel> t) {
        super(t);
    }

    @Override
    public void serialize(SimulationModel value, JsonGenerator gen, SerializerProvider provider) throws IOException {
        gen.writeStartObject();

        // Write "classes" field
        Map<Class, List<InstanceProperty>> classes = value.getClasses();
        gen.writeFieldName("classes");
        gen.writeStartObject();
        for (Map.Entry<Class, List<InstanceProperty>> entry : classes.entrySet()) {
            gen.writeFieldName(entry.getKey().getSimpleName());

            gen.writeStartObject();
            for(InstanceProperty property : entry.getValue()) {
                gen.writeObjectField(property.getName(), property);
            }
            gen.writeEndObject();

        }
        gen.writeEndObject();

        // Write "tree" field
        gen.writeObjectField("tree", value.getTree());

        gen.writeEndObject();


    }
}
