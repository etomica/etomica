package etomica.server;

import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.databind.module.SimpleModule;
import com.google.inject.AbstractModule;
import com.google.inject.Injector;
import com.google.inject.Provides;
import etomica.data.IEtomicaDataSource;
import etomica.meta.ComponentIndex;
import etomica.meta.DataSourceIndex;
import etomica.meta.SimulationModel;
import etomica.meta.properties.Property;
import etomica.meta.wrappers.Wrapper;
import etomica.server.health.BasicHealthCheck;
import etomica.server.resources.*;
import etomica.server.resources.data.DataStreamWebsocket;
import etomica.server.serializers.DataSerializer;
import etomica.server.serializers.PropertySerializer;
import etomica.server.serializers.SimulationModelSerializer;
import etomica.server.serializers.WrapperSerializer;
import etomica.simulation.Simulation;
import io.dropwizard.Application;
import io.dropwizard.jackson.Jackson;
import io.dropwizard.setup.Bootstrap;
import io.dropwizard.setup.Environment;
import io.dropwizard.websockets.WebsocketBundle;
import org.eclipse.jetty.servlets.CrossOriginFilter;
import org.eclipse.jetty.websocket.jsr356.server.BasicServerEndpointConfig;
import ru.vyarus.dropwizard.guice.GuiceBundle;

import javax.inject.Inject;
import javax.inject.Singleton;
import javax.servlet.DispatcherType;
import javax.servlet.FilterRegistration;
import javax.websocket.server.ServerEndpointConfig;
import java.util.EnumSet;
import java.util.Map;
import java.util.Timer;
import java.util.UUID;
import java.util.concurrent.ConcurrentHashMap;

public class EtomicaServer extends Application<EtomicaServerConfig> {
    private final Map<UUID, SimulationModel> simStore = new ConcurrentHashMap<>();
    private final Timer timer = new Timer();
    private final ObjectMapper mapper = Jackson.newObjectMapper();

    {
        SimpleModule mod = new SimpleModule("Etomica Module");
        mod.addSerializer(new PropertySerializer(Property.class));
        mod.addSerializer(new WrapperSerializer(Wrapper.class));
        mod.addSerializer(new SimulationModelSerializer(SimulationModel.class));
        mod.addSerializer(new DataSerializer());
        mapper.registerModule(mod);
    }

    @Override
    public String getName() {
        return "EtomicaServer";
    }

    @Override
    public void initialize(Bootstrap<EtomicaServerConfig> bootstrap) {
        bootstrap.setObjectMapper(mapper);

        bootstrap.addBundle(GuiceBundle.builder()
                .enableAutoConfig(getClass().getPackage().getName())
                .modules(new WebSocketModule(), new EtomicaServerModule(bootstrap.getObjectMapper()))
                .build()
        );

        WebsocketBundle wsBundle = new WebsocketBundle(new WSConfigurator());
        wsBundle.addEndpoint(EchoServer.class);
        wsBundle.addEndpoint(ConfigurationStreamResource.class);
        wsBundle.addEndpoint(DataStreamWebsocket.class);

        bootstrap.addBundle(wsBundle);

        super.initialize(bootstrap);
    }

    @Override
    public void run(EtomicaServerConfig configuration, Environment environment) throws Exception {
        environment.healthChecks().register("basic", new BasicHealthCheck());
        configureCors(environment);

//        environment.jersey().register(new SimulationsIndexResource(new ComponentIndex<>(Simulation.class)));
//        environment.jersey().register(new SimulationResource(simStore));
//        environment.jersey().register(new ControlResource(simStore));

    }

    public static void main(String[] args) throws Exception {
        new EtomicaServer().run(args);
    }

    private static void configureCors(Environment env) {
        final FilterRegistration.Dynamic cors = env.servlets().addFilter("CORS", CrossOriginFilter.class);

        cors.setInitParameter(CrossOriginFilter.ALLOWED_ORIGINS_PARAM, "*");
        cors.setInitParameter(CrossOriginFilter.ALLOWED_HEADERS_PARAM, "X-Requested-With,Content-Type,Accept,Origin,Authorization");
        cors.setInitParameter(CrossOriginFilter.ALLOWED_METHODS_PARAM, "OPTIONS,GET,PUT,POST,DELETE,HEAD");
        cors.setInitParameter(CrossOriginFilter.ALLOW_CREDENTIALS_PARAM, "true");

        // Add URL mapping
        cors.addMappingForUrlPatterns(EnumSet.allOf(DispatcherType.class), true, "/*");
    }

    private static class WSConfigurator extends ServerEndpointConfig.Configurator {
        @Inject
        private static Injector injector;

        @Override
        public <T> T getEndpointInstance(Class<T> endpointClass) throws InstantiationException {
            return injector.getInstance(endpointClass);
        }
    }

    private static class WebSocketModule extends AbstractModule {

        @Override
        protected void configure() {
            requestStaticInjection(WSConfigurator.class);
        }
    }

    private static class EtomicaServerModule extends AbstractModule {

        private final ObjectMapper mapper;

        public EtomicaServerModule(ObjectMapper mapper) {
            this.mapper = mapper;
        }

        @Override
        protected void configure() {
        }

        @Provides @Singleton
        ComponentIndex<Simulation> provideSimulationIndex() {
            return new ComponentIndex<>(Simulation.class);
        }

        @Provides @Singleton
        DataSourceIndex provideDataSourceIndex() {
            return new DataSourceIndex(IEtomicaDataSource.class);
        }

        @Provides @Singleton
        Timer provideTimer() {
            return new Timer();
        }

        @Provides @Singleton
        ObjectMapper provideObjectMapper() {
            return mapper;
        }
    }
}
