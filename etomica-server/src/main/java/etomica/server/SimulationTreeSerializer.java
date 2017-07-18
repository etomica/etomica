package etomica.server;

import com.fasterxml.jackson.core.JsonGenerator;
import com.fasterxml.jackson.databind.SerializerProvider;
import com.fasterxml.jackson.databind.ser.std.StdSerializer;
import etomica.meta.InstanceProperty;

import java.io.IOException;
import java.util.Map;
import java.util.Map.Entry;

public class SimulationTreeSerializer extends StdSerializer<SimulationTree> {

    public SimulationTreeSerializer(Class<SimulationTree> t) {
        super(t);
    }

    @Override
    public void serialize(SimulationTree value, JsonGenerator gen, SerializerProvider provider) throws IOException {
        gen.writeStartObject();
        gen.writeStringField("class", value.getData().getWrappedClass().getName());
        for(Entry<String, Object> e : value.getData().getValues().entrySet()) {
            gen.writeObjectField(e.getKey(), e.getValue());
        }
        for(SimulationTree tree : value.getChildren()) {
            if(tree.getData() == null) {
                gen.writeNullField(tree.getPropName());
            } else {
                gen.writeObjectField(tree.getPropName(), tree);
            }
        }
        gen.writeEndObject();
    }


}
