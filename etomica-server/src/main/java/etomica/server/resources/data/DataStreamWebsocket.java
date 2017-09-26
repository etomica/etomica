package etomica.server.resources.data;

import com.codahale.metrics.annotation.Metered;
import com.codahale.metrics.annotation.Timed;
import etomica.data.DataDump;
import etomica.data.DataPumpListener;
import etomica.data.IData;
import etomica.data.IDataSource;
import etomica.meta.SimulationModel;
import etomica.server.dao.DataStreamStore;
import etomica.server.dao.SimulationStore;
import etomica.simulation.Simulation;

import javax.inject.Inject;
import javax.websocket.OnOpen;
import javax.websocket.Session;
import javax.websocket.server.PathParam;
import javax.websocket.server.ServerEndpoint;
import java.util.Timer;
import java.util.TimerTask;
import java.util.UUID;

@ServerEndpoint(
        value="/simulations/{simId}/data/{dataId}"
)
@Metered
@Timed
public class DataStreamWebsocket {
    private final SimulationStore simStore;
    private final DataStreamStore dataStore;
    private final Timer timer;

    @Inject
    public DataStreamWebsocket(SimulationStore simStore, DataStreamStore dataStore, Timer timer) {
        this.simStore = simStore;
        this.dataStore = dataStore;
        this.timer = timer;
    }

    @OnOpen
    public void onOpen(final Session session, @PathParam("simId") String simId, @PathParam("dataId") String dataId) {
        SimulationModel model = simStore.get(UUID.fromString(simId));
        DataStreamStore.DataPlumbing dataPlumbing = dataStore.get(UUID.fromString(dataId));

        model.getSimulation().getIntegrator().getEventManager().addListener(dataPlumbing.getPump());
        timer.schedule(new DataStreamTimerTask(session, dataPlumbing.getDump(), model.getSimulation()), 0, 333);
    }

    public static class DataStreamTimerTask extends TimerTask {
        private final Session session;
        private final IDataSource dump;
        private final Simulation sim;

        private DataStreamTimerTask(Session session, IDataSource dump, Simulation sim) {
            this.session = session;
            this.dump = dump;
            this.sim = sim;
        }

        @Override
        public void run() {
            if(sim.getController().isPaused() || !sim.getController().isActive()) {
                return;
            }

            sim.getController().doActionNow(() -> {
                IData data = dump.getData();
                if(data != null) {
                    session.getAsyncRemote().sendObject(data.getValue(0));
                }
            });
        }
    }
}
