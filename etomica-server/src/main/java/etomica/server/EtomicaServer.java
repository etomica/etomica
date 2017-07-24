package etomica.server;

import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.databind.module.SimpleModule;
import etomica.meta.properties.Property;
import etomica.meta.SimulationModel;
import etomica.meta.wrappers.Wrapper;
import etomica.simulation.prototypes.HSMD2D;
import org.jooby.Jooby;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;

public class EtomicaServer extends Jooby {
    public static final ObjectMapper mapper = new ObjectMapper();
    private static final SimpleModule mod = new SimpleModule("Etomica Module");
    static {
        mod.addSerializer(new PropertySerializer(Property.class));
        mod.addSerializer(new WrapperSerializer(Wrapper.class));
        mod.addSerializer(new SimulationModelSerializer(SimulationModel.class));
        mapper.registerModule(mod);
    }


    {
        get("/", () -> "Hello World");
    }

    public static void main(String[] args) throws IOException {
//        ComponentIndex<Simulation> idx = new ComponentIndex<>(Simulation.class);
//        Set<Class<? extends Simulation>> set = idx.getComponentSet();

//        SimulationWrapper simWrapper = new SimulationWrapper(new HSMD2D());
//        SimulationModel model = new SimulationModel(simWrapper);
        SimulationModel model = new SimulationModel(new HSMD2D());
        String s = mapper.writerWithDefaultPrettyPrinter().writeValueAsString(model);
        Files.write(Paths.get("./simStructure.json"), s.getBytes());

//        run(EtomicaServer::new, args);
    }
}
