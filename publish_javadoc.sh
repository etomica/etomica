#!/usr/bin/env bash

if [ ! "$TRAVIS" == "true" ]; then
    echo -en "\e[31m"
    echo -n "This script should only be executed by the Travis build server! Exiting."
    echo -e "\e[39m"
    exit 1
fi

if [ "$TRAVIS_REPO_SLUG" == "etomica/etomica" ] && [ "$TRAVIS_JDK_VERSION" == "oraclejdk8" ] && [ "$TRAVIS_PULL_REQUEST" == "false" ] && [ "$TRAVIS_BRANCH" == "master" ]; then
    echo "Generating javadoc..."
    ./gradlew javadocAll || { echo "Javadoc build failed!"; exit 1; }

    echo "Publishing javadoc..."
    cd $HOME
    git config --global user.email "travis@travis-ci.org"
    git config --global user.name "travis-ci"
    git clone https://${GH_TOKEN}@github.com/etomica/javadoc etomica-javadoc

    cd etomica-javadoc
    rsync -r --delete --exclude=.git ${TRAVIS_BUILD_DIR}/docs/javadoc/ .

    [ -z "$(git status --porcelain)" ] && { echo "No changes to javadoc"; exit 0; }

    COMMIT_MSG="Latest javadoc auto-pushed by travis

    travis build ${TRAVIS_BUILD_NUMBER}
    from commit ${TRAVIS_REPO_SLUG}@${TRAVIS_COMMIT}"

    git add .
    git commit -m "$COMMIT_MSG"
    git push origin master
    echo "Successfully pushed javadoc to https://github.com/etomica/javadoc"

else
    echo "Not on master branch of main repository, exiting."
fi
