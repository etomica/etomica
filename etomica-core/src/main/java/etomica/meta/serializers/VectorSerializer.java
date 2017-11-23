package etomica.meta.serializers;

import com.fasterxml.jackson.core.JsonGenerator;
import com.fasterxml.jackson.databind.SerializerProvider;
import com.fasterxml.jackson.databind.ser.std.StdSerializer;
import etomica.space.Vector;

import java.io.IOException;

public class VectorSerializer extends StdSerializer<Vector> {
    protected VectorSerializer() {
        super(Vector.class);
    }

    @Override
    public void serialize(Vector value, JsonGenerator gen, SerializerProvider provider) throws IOException {
        double[] arr = value.toArray();
        gen.writeArray(arr, 0, arr.length);

    }
}
