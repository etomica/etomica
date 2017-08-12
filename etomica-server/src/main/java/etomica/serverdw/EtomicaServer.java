package etomica.serverdw;

import com.fasterxml.jackson.databind.module.SimpleModule;
import etomica.meta.ComponentIndex;
import etomica.meta.SimulationModel;
import etomica.meta.properties.Property;
import etomica.meta.wrappers.Wrapper;
import etomica.serverdw.health.BasicHealthCheck;
import etomica.serverdw.resources.EchoServer;
import etomica.serverdw.resources.SimulationResource;
import etomica.serverdw.resources.SimulationsIndexResource;
import etomica.serverdw.serializers.PropertySerializer;
import etomica.serverdw.serializers.SimulationModelSerializer;
import etomica.serverdw.serializers.WrapperSerializer;
import etomica.simulation.Simulation;
import io.dropwizard.Application;
import io.dropwizard.setup.Bootstrap;
import io.dropwizard.setup.Environment;
import io.dropwizard.websockets.WebsocketBundle;

import java.util.HashMap;
import java.util.concurrent.ConcurrentHashMap;

public class EtomicaServer extends Application<EtomicaServerConfig> {
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
        environment.jersey().register(new SimulationResource(new ConcurrentHashMap<>()));

    }

    public static void main(String[] args) throws Exception {
        new EtomicaServer().run(args);
    }
}
