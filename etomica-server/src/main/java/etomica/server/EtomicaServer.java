package etomica.server;

import com.fasterxml.jackson.databind.module.SimpleModule;
import etomica.meta.ComponentIndex;
import etomica.meta.SimulationModel;
import etomica.meta.properties.Property;
import etomica.meta.wrappers.Wrapper;
import etomica.server.health.BasicHealthCheck;
import etomica.server.resources.ControlResource;
import etomica.server.resources.EchoServer;
import etomica.server.resources.SimulationResource;
import etomica.server.resources.SimulationsIndexResource;
import etomica.server.serializers.PropertySerializer;
import etomica.server.serializers.SimulationModelSerializer;
import etomica.server.serializers.WrapperSerializer;
import etomica.simulation.Simulation;
import io.dropwizard.Application;
import io.dropwizard.setup.Bootstrap;
import io.dropwizard.setup.Environment;
import io.dropwizard.websockets.WebsocketBundle;

import java.util.Map;
import java.util.UUID;
import java.util.concurrent.ConcurrentHashMap;

public class EtomicaServer extends Application<EtomicaServerConfig> {
    private final Map<UUID, SimulationModel> simStore = new ConcurrentHashMap<>();

    @Override
    public String getName() {
        return "EtomicaServer";
    }

    @Override
    public void initialize(Bootstrap<EtomicaServerConfig> bootstrap) {

        SimpleModule mod = new SimpleModule("Etomica Module");
        mod.addSerializer(new PropertySerializer(Property.class));
        mod.addSerializer(new WrapperSerializer(Wrapper.class));
        mod.addSerializer(new SimulationModelSerializer(SimulationModel.class));

        bootstrap.getObjectMapper().registerModule(mod);

        bootstrap.addBundle(new WebsocketBundle(EchoServer.class));

        super.initialize(bootstrap);
    }

    @Override
    public void run(EtomicaServerConfig configuration, Environment environment) throws Exception {
        environment.healthChecks().register("basic", new BasicHealthCheck());

        environment.jersey().register(new SimulationsIndexResource(new ComponentIndex<>(Simulation.class)));
        environment.jersey().register(new SimulationResource(simStore));
        environment.jersey().register(new ControlResource(simStore));

    }

    public static void main(String[] args) throws Exception {
        new EtomicaServer().run(args);
    }
}
