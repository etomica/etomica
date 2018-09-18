#!/usr/bin/env bash

# fail on any non-zero exit codes
set -e

if [ ! "$TRAVIS" == "true" ]; then
    echo "$(tput setaf 1)This script should only be executed by the Travis build server! Exiting.$(tput sgr 0)"
    exit 1
fi

# check for: etomica repo (not a fork), consistent jdk version (in case we ever test against multiple),
#            not a pull request, and only commits to master branch.
if [ "$TRAVIS_REPO_SLUG" == "etomica/etomica" ] && [ "$TRAVIS_JDK_VERSION" == "oraclejdk9" ] && [ "$TRAVIS_PULL_REQUEST" == "false" ] && [ "$TRAVIS_BRANCH" == "master" ]; then
    COMMIT_AUTHOR_EMAIL="$(git --no-pager show -s --format='%ce' ${TRAVIS_COMMIT})"
    echo "Generating javadoc..."
    ./gradlew javadocAll || { echo "Javadoc build failed!"; exit 1; }


    echo "Publishing javadoc..."
    cd $HOME
    git config --global user.email "${COMMIT_AUTHOR_EMAIL}"
    git config --global user.name "travis-ci"
    git clone https://${GH_TOKEN}@github.com/etomica/javadoc etomica-javadoc

    cd etomica-javadoc
    rsync -rh --stats --delete --exclude=.git ${TRAVIS_BUILD_DIR}/docs/javadoc/ .

    [ -z "$(git status --porcelain)" ] && { echo "No changes to javadoc"; exit 0; }

    COMMIT_MSG="Latest javadoc auto-pushed by travis

    travis build ${TRAVIS_BUILD_NUMBER}
    from commit ${TRAVIS_REPO_SLUG}@${TRAVIS_COMMIT}"

    git add -A
    git commit -m "$COMMIT_MSG"
    echo "Pushing to javadoc repo..."
    git push --quiet origin master > /dev/null 2>&1
    echo "Successfully pushed javadoc to https://github.com/etomica/javadoc"

else
    echo "Not on master branch of main repository, exiting."
fi
