package etomica.server;

import com.fasterxml.jackson.databind.module.SimpleModule;
import etomica.meta.ComponentIndex;
import etomica.meta.SimulationModel;
import etomica.meta.properties.Property;
import etomica.meta.wrappers.Wrapper;
import etomica.server.health.BasicHealthCheck;
import etomica.server.resources.*;
import etomica.server.serializers.PropertySerializer;
import etomica.server.serializers.SimulationModelSerializer;
import etomica.server.serializers.WrapperSerializer;
import etomica.simulation.Simulation;
import io.dropwizard.Application;
import io.dropwizard.setup.Bootstrap;
import io.dropwizard.setup.Environment;
import io.dropwizard.websockets.WebsocketBundle;

import javax.websocket.server.ServerEndpointConfig;
import java.util.Map;
import java.util.Timer;
import java.util.UUID;
import java.util.concurrent.ConcurrentHashMap;

public class EtomicaServer extends Application<EtomicaServerConfig> {
    private final Map<UUID, SimulationModel> simStore = new ConcurrentHashMap<>();
    private final Timer timer = new Timer();

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

        final ServerEndpointConfig wsConfig = ServerEndpointConfig.Builder
                .create(ConfigurationStreamResource.class, "/simulations/{id}/configuration")
                .build();

        wsConfig.getUserProperties().put("simStore", simStore);
        wsConfig.getUserProperties().put("timer", timer);

        WebsocketBundle wsBundle = new WebsocketBundle(wsConfig);
        wsBundle.addEndpoint(EchoServer.class);

        bootstrap.addBundle(wsBundle);

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
