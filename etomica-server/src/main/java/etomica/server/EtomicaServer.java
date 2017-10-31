package etomica.server;

import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.databind.ObjectWriter;
import com.fasterxml.jackson.databind.module.SimpleModule;
import com.google.inject.AbstractModule;
import com.google.inject.Injector;
import com.google.inject.Provides;
import etomica.data.IDataSource;
import etomica.meta.ComponentIndex;
import etomica.meta.DataSourceIndex;
import etomica.server.health.BasicHealthCheck;
import etomica.server.resources.ConfigurationWebsocket;
import etomica.server.resources.EchoServer;
import etomica.server.resources.data.DataStreamWebsocket;
import etomica.server.serializers.DataSerializer;
import etomica.server.serializers.WrapperSerializer;
import etomica.simulation.Simulation;
import io.dropwizard.Application;
import io.dropwizard.jackson.Jackson;
import io.dropwizard.setup.Bootstrap;
import io.dropwizard.setup.Environment;
import io.dropwizard.websockets.WebsocketBundle;
import org.eclipse.jetty.servlets.CrossOriginFilter;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import ru.vyarus.dropwizard.guice.GuiceBundle;

import javax.inject.Inject;
import javax.inject.Singleton;
import javax.servlet.DispatcherType;
import javax.servlet.FilterRegistration;
import javax.websocket.Extension;
import javax.websocket.HandshakeResponse;
import javax.websocket.server.HandshakeRequest;
import javax.websocket.server.ServerEndpointConfig;
import java.util.*;

public class EtomicaServer extends Application<EtomicaServerConfig> {
    private static final ObjectMapper mapper = Jackson.newObjectMapper();
    private static final Logger log = LoggerFactory.getLogger(EtomicaServer.class);

    static {
        SimpleModule mod = new SimpleModule("Etomica Module");
        mod.addSerializer(new WrapperSerializer());
        mod.addSerializer(new DataSerializer());
        mapper.registerModule(mod);
    }

    public static ObjectWriter objectWriter() {
        return mapper.writer();
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

        WSConfigurator wsConfigurator = new WSConfigurator();
        WebsocketBundle wsBundle = new WebsocketBundle(wsConfigurator);
        wsBundle.addEndpoint(EchoServer.class);
        wsBundle.addEndpoint(ConfigurationWebsocket.class);
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
        public List<Extension> getNegotiatedExtensions(List<Extension> installed, List<Extension> requested) {
            // permessage-deflate is causing problems with the deflater being closed. Nuke it until I can figure out how to fix it.
            return Collections.emptyList();
        }

        @Override
        public void modifyHandshake(ServerEndpointConfig sec, HandshakeRequest request, HandshakeResponse response) {
        }

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
            return new DataSourceIndex(IDataSource.class);
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
