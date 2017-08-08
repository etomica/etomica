package etomica.server;

import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.databind.module.SimpleModule;
import etomica.meta.ComponentIndex;
import etomica.meta.properties.Property;
import etomica.meta.SimulationModel;
import etomica.meta.wrappers.Wrapper;
import etomica.server.controllers.SimulationController;
import etomica.server.controllers.SimulationIndexController;
import etomica.simulation.Simulation;
import etomica.simulation.prototypes.HSMD2D;
import org.jooby.Jooby;
import org.jooby.MediaType;
import org.jooby.handlers.Cors;
import org.jooby.handlers.CorsHandler;
import org.jooby.json.Jackson;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;

public class EtomicaServer extends Jooby {
    private static final ComponentIndex<Simulation> simIndex = new ComponentIndex<>(Simulation.class);

    public static final ObjectMapper mapper = new ObjectMapper();
    private static final SimpleModule mod = new SimpleModule("Etomica Module");
    static {
        mod.addSerializer(new PropertySerializer(Property.class));
        mod.addSerializer(new WrapperSerializer(Wrapper.class));
        mod.addSerializer(new SimulationModelSerializer(SimulationModel.class));
        mapper.registerModule(mod);
    }

    public static final Map<UUID, SimulationModel> simulations = new ConcurrentHashMap<>();




    {
        use(new Jackson(mapper));
        use("*", new CorsHandler(new Cors()));
        use(SimulationIndexController.class);
        use(SimulationController.class);

//        get("/simulation/:id", (req, resp) -> {
//
//            UUID id = UUID.fromString(req.param("id").value());
//            SimulationModel model = simulations.get(id);
//            resp.type(MediaType.json).send(model);
//        });
    }

    public static void main(String[] args) throws IOException {
//        ComponentIndex<Simulation> idx = new ComponentIndex<>(Simulation.class);
//        Set<Class<? extends Simulation>> set = idx.getComponentSet();

//        SimulationWrapper simWrapper = new SimulationWrapper(new HSMD2D());
//        SimulationModel model = new SimulationModel(simWrapper);
//        SimulationModel model = new SimulationModel(new HSMD2D());
//        String s = mapper.writerWithDefaultPrettyPrinter().writeValueAsString(model);
//        Files.write(Paths.get("./simStructure.json"), s.getBytes());

        run(EtomicaServer::new, args);
    }
}
