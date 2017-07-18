package etomica.server;

import com.fasterxml.jackson.core.JsonProcessingException;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.databind.module.SimpleModule;
import etomica.meta.ComponentIndex;
import etomica.meta.InstanceProperty;
import etomica.meta.wrappers.SimulationWrapper;
import etomica.meta.wrappers.Wrapper;
import etomica.simulation.Simulation;
import etomica.simulation.prototypes.HSMD2D;
import org.jooby.Jooby;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Set;

public class EtomicaServer extends Jooby {
    public static final ObjectMapper mapper = new ObjectMapper();
    private static final SimpleModule mod = new SimpleModule("Etomica Module");
    static {
        mod.addSerializer(new PropertySerializer(InstanceProperty.class));
        mod.addSerializer(new WrapperSerializer(Wrapper.class));
        mod.addSerializer(new SimulationTreeSerializer(SimulationTree.class));
        mod.addSerializer(new SimulationModelSerializer(SimulationModel.class));
        mapper.registerModule(mod);
    }


    {
        get("/", () -> "Hello World");
    }

    public static void main(String[] args) throws IOException {
//        ComponentIndex<Simulation> idx = new ComponentIndex<>(Simulation.class);
//        Set<Class<? extends Simulation>> set = idx.getComponentSet();

        SimulationWrapper simWrapper = new SimulationWrapper(new HSMD2D());
        SimulationModel model = new SimulationModel(simWrapper);
        String s = mapper.writerWithDefaultPrettyPrinter().writeValueAsString(model);
        Files.write(Paths.get("./simStructure.json"), s.getBytes());

//        run(EtomicaServer::new, args);
    }
}
