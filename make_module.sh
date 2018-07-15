#!/usr/bin/env bash

./gradlew :etomica-modules:jar && jlink --module-path $PWD/etomica-modules/build/libs/etomica-modules.jar:$PWD/etomica-apps/build/libs/etomica-apps.jar:$PWD/etomica-core/build/libs/etomica-core.jar:$PWD/libs/ptplot5.10.jar:$PWD/vendor/Jama-1.0.3.jar:$PWD/vendor/commons-math3-3.3.jar:$PWD/etomica-graph/build/libs/etomica-graph.jar:$PWD/etomica-graphics3D/build/libs/etomica-graphics3D.jar:$JAVA_HOME/jmods \
    --add-modules etomica.modules --launcher stats=etomica.modules/etomica.modules.statistics.StatisticsMCGraphic --output stats_module --strip-debug --compress 2 --no-header-files --no-man-pages
