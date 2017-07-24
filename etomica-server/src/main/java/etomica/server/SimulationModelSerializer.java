package etomica.server;

import com.fasterxml.jackson.core.JsonGenerator;
import com.fasterxml.jackson.databind.SerializerProvider;
import com.fasterxml.jackson.databind.ser.std.StdSerializer;
import etomica.meta.SimulationModel;
import etomica.meta.wrappers.Wrapper;

import java.io.IOException;
import java.util.Collection;

public class SimulationModelSerializer extends StdSerializer<SimulationModel> {


    protected SimulationModelSerializer(Class<SimulationModel> t) {
        super(t);
    }

    @Override
    public void serialize(SimulationModel value, JsonGenerator gen, SerializerProvider provider) throws IOException {
//        gen.writeStartObject();

        // Write "classes" field
//        Map<Class, List<Property>> classes = value.getClasses();
//        gen.writeFieldName("classes");

        gen.writeObject(value.allWrappers());

//        Collection<Wrapper> allWrappers = value.allWrappers();
//        gen.writeStartObject();
//        for (Wrapper wrapper : allWrappers) {
//            gen.writeFieldName(entry.getKey().getSimpleName());

//            gen.writeStartObject();
//            gen.writeObject(wrapper);
//            for(Property property : entry.getValue()) {
//                gen.writeObjectField(property.getName(), property);
//            }
//            gen.writeEndObject();

//        }
//        gen.writeEndObject();

//        gen.writeEndObject();


    }
}
