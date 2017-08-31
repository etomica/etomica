package etomica.server.serializers;

import com.fasterxml.jackson.core.JsonGenerator;
import com.fasterxml.jackson.databind.SerializerProvider;
import com.fasterxml.jackson.databind.ser.std.StdSerializer;
import etomica.meta.wrappers.Wrapper;
import etomica.space.Vector;

import java.io.IOException;

public class WrapperSerializer extends StdSerializer<Wrapper<?>> {

    public WrapperSerializer(Class<Wrapper> t) {
        super(t, true);
    }

    @Override
    public void serialize(Wrapper<?> value, JsonGenerator gen, SerializerProvider provider) throws IOException {
        // TODO: maybe find a cleaner way to do this
        if(Vector.class.isAssignableFrom(value.getWrappedClass())) {
            // value is a vector, just write the coordinates
            gen.writeObject(((Vector) value.getWrapped()).toArray());
            return;
        }
        gen.writeStartObject();
        gen.writeStringField("class", value.getWrappedClass().getSimpleName());
        gen.writeNumberField("id", value.getWrappedId());
        value.getValues().forEach((name, obj) -> {
            try {
                gen.writeObjectField(name, obj);
            } catch (IOException e) {
                e.printStackTrace();
            }
        });
//        for(Property p : value.getProperties()) {
//            gen.writeObjectField(p.getName(), p);
//        }
        gen.writeEndObject();
    }
}
