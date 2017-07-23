package etomica.server;

import com.fasterxml.jackson.core.JsonGenerator;
import com.fasterxml.jackson.databind.SerializerProvider;
import com.fasterxml.jackson.databind.ser.std.StdSerializer;
import etomica.meta.wrappers.CollectionWrapper;

import java.io.IOException;
import java.util.Map.Entry;

public class SimulationTreeSerializer extends StdSerializer<SimulationTree> {

    public SimulationTreeSerializer(Class<SimulationTree> t) {
        super(t);
    }

    @Override
    public void serialize(SimulationTree value, JsonGenerator gen, SerializerProvider provider) throws IOException {

        if(value.getWrapper() instanceof CollectionWrapper<?>) {
            gen.writeObject(value.getChildren());
        } else {
            gen.writeStartObject();
            gen.writeStringField("class", value.getWrapper().getWrappedClass().getName());
            for(Entry<String, Object> e : value.getWrapper().getValues().entrySet()) {
                gen.writeObjectField(e.getKey(), e.getValue());
            }
            for(SimulationTree tree : value.getChildren()) {
                if(tree.getWrapper() == null) {
                    gen.writeNullField(((PropertyTree) tree).getPropName());
                } else {
                    gen.writeObjectField(((PropertyTree) tree).getPropName(), tree);
                }
            }
            gen.writeEndObject();
        }

    }


}
