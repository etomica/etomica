package etomica.server.resources.data;

import com.codahale.metrics.annotation.Metered;
import com.codahale.metrics.annotation.Timed;
import com.fasterxml.jackson.core.JsonProcessingException;
import com.fasterxml.jackson.databind.ObjectMapper;
import etomica.data.*;
import etomica.meta.SimulationModel;
import etomica.server.dao.DataStreamStore;
import etomica.server.dao.SimulationStore;
import etomica.server.representations.DataAndInfo;
import etomica.simulation.Simulation;

import javax.inject.Inject;
import javax.websocket.*;
import javax.websocket.server.PathParam;
import javax.websocket.server.ServerEndpoint;
import java.util.Timer;
import java.util.TimerTask;
import java.util.UUID;

@ServerEndpoint(
        value="/simulations/{simId}/data/{dataId}",
        encoders = {DataStreamWebsocket.DataAndInfoEncoder.class}
)
@Metered
@Timed
public class DataStreamWebsocket {
    private final SimulationStore simStore;
    private final DataStreamStore dataStore;
    private final Timer timer;
    private final ObjectMapper mapper;

    @Inject
    public DataStreamWebsocket(SimulationStore simStore, DataStreamStore dataStore, Timer timer, ObjectMapper mapper) {
        this.simStore = simStore;
        this.dataStore = dataStore;
        this.timer = timer;
        this.mapper = mapper;
    }

    @OnOpen
    public void onOpen(final Session session, @PathParam("simId") String simId, @PathParam("dataId") String dataId) {
        session.getUserProperties().put("mapper", mapper);
        SimulationModel model = simStore.get(UUID.fromString(simId));
        DataStreamStore.DataPlumbing dataPlumbing = dataStore.get(UUID.fromString(dataId));

        model.getSimulation().getIntegrator().getEventManager().addListener(dataPlumbing.getPump());
        timer.schedule(new DataStreamTimerTask(session, dataPlumbing.getDump(), model.getSimulation()), 0, 333);
    }

    public static class DataStreamTimerTask extends TimerTask {
        private final Session session;
        private final IDataSource dump;
        private final Simulation sim;
        private final DataAndInfo dataAndInfo;

        private DataStreamTimerTask(Session session, IDataSource dump, Simulation sim) {
            this.session = session;
            this.dump = dump;
            this.sim = sim;
            this.dataAndInfo = new DataAndInfo();
        }

        @Override
        public void run() {
            if(sim.getController().isPaused() || !sim.getController().isActive()) {
                return;
            }

            sim.getController().doActionNow(() -> {
                IData data = dump.getData();
                this.dataAndInfo.setData(dump.getData());
                this.dataAndInfo.setDataInfo(dump.getDataInfo());
                if(data != null) {
                     session.getAsyncRemote().sendObject(this.dataAndInfo);
                }
            });
        }
    }

    public static class DataAndInfoEncoder implements Encoder.Text<DataAndInfo> {
        private EndpointConfig config;

        @Override
        public String encode(DataAndInfo object) throws EncodeException {
            ObjectMapper mapper = (ObjectMapper) config.getUserProperties().get("mapper");
            try {
                return mapper.writeValueAsString(object);
            } catch (JsonProcessingException e) {
                throw new EncodeException(object, "JSON error", e);
            }
        }

        @Override
        public void init(EndpointConfig config) {
            this.config = config;
        }

        @Override
        public void destroy() {

        }
    }

}
