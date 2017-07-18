package etomica.server;

import com.fasterxml.jackson.core.JsonGenerator;
import com.fasterxml.jackson.databind.SerializerProvider;
import com.fasterxml.jackson.databind.ser.std.StdSerializer;
import etomica.meta.InstanceProperty;
import etomica.meta.wrappers.Wrapper;

import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class WrapperSerializer extends StdSerializer<Wrapper<?>> {

    protected WrapperSerializer(Class<Wrapper> t) {
        super(t, true);
    }

    @Override
    public void serialize(Wrapper<?> value, JsonGenerator gen, SerializerProvider provider) throws IOException {
        gen.writeStartObject();
        gen.writeStringField("class", value.getWrappedClass().getSimpleName());
        for(InstanceProperty p : value.getProperties()) {
            gen.writeObjectField(p.getName(), p);
        }
        gen.writeEndObject();
    }
}
