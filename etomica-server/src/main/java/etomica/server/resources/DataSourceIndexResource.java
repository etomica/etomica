package etomica.server.resources;

import com.google.common.base.Predicate;
import com.google.common.base.Predicates;
import etomica.box.Box;
import etomica.data.IDataSink;
import etomica.data.IEtomicaDataSource;
import etomica.meta.DataSourceIndex;
import org.reflections.ReflectionUtils;

import javax.inject.Inject;
import javax.ws.rs.GET;
import javax.ws.rs.Path;
import javax.ws.rs.Produces;
import javax.ws.rs.core.MediaType;
import java.lang.reflect.Constructor;
import java.util.List;
import java.util.stream.Collectors;

import static org.reflections.ReflectionUtils.withParametersAssignableTo;
import static org.reflections.ReflectionUtils.withParametersCount;
import static org.reflections.ReflectionUtils.withTypeAssignableTo;

@Path("/simulations/{simId}/data")
@Produces(MediaType.APPLICATION_JSON)
public class DataSourceIndexResource {
    private final DataSourceIndex dataSourceIndex;

    @Inject
    public DataSourceIndexResource(DataSourceIndex index) {
        this.dataSourceIndex = index;
    }

    @GET
    public List<String> list() {
        return dataSourceIndex.getComponentSet().stream()
                .filter(cls -> !IDataSink.class.isAssignableFrom(cls))
                .filter(cls -> {
                    return ReflectionUtils.getAllConstructors(cls, Predicates.or(
                            withParametersCount(0),
                            withParametersAssignableTo(Box.class)
                    )).size() != 0;
                })
                .map(Class::getCanonicalName)
                .collect(Collectors.toList());
    }

}
