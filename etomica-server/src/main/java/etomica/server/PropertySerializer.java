package etomica.server;

import com.fasterxml.jackson.core.JsonGenerator;
import com.fasterxml.jackson.databind.SerializerProvider;
import com.fasterxml.jackson.databind.ser.std.StdSerializer;
import etomica.meta.InstanceProperty;

import java.io.IOException;

public class PropertySerializer extends StdSerializer<InstanceProperty>{

    protected PropertySerializer(Class<InstanceProperty> t) {
        super(t);
    }

    @Override
    public void serialize(InstanceProperty value, JsonGenerator gen, SerializerProvider provider) throws IOException {

        gen.writeStartObject();

        gen.writeStringField("type", value.getPropertyType().getSimpleName());
        gen.writeBooleanField("indexed", value.isIndexedProperty());

        gen.writeArrayFieldStart("methods");
        if(value.canRead()) { gen.writeString("read"); }
        if(value.canWrite()) { gen.writeString("write"); }
        if(value.canAdd()) { gen.writeString("add"); }
        if(value.canRemove()) { gen.writeString("remove"); }
        if(value.canCount()) { gen.writeString("count"); }
        gen.writeEndArray();

        gen.writeEndObject();

    }
}
