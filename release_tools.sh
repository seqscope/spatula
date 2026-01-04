#!/bin/bash

# 1. SAFETY CHECK: Ensure we are on the 'main' branch
CURRENT_BRANCH=$(git rev-parse --abbrev-ref HEAD)
if [ "$CURRENT_BRANCH" != "main" ]; then
  echo "Error: Releases must be performed from the 'main' branch."
  echo "Current branch: $CURRENT_BRANCH"
  exit 1
fi

# 2. SAFETY CHECK: Ensure the working directory is clean
# (You don't want to release code you haven't committed yet)
if [ -n "$(git status --porcelain)" ]; then
  echo "Error: Working directory is not clean. Commit or stash changes first."
  exit 1
fi

# 1. EXTRACT VERSION from CMakeLists.txt
# This looks for the line "project(... VERSION X.Y.Z ...)"
VERSION=$(grep -oP 'project.*VERSION \K[0-9]+\.[0-9]+\.[0-9]+' CMakeLists.txt)

if [ -z "$VERSION" ]; then
    echo "Error: Could not find version in CMakeLists.txt"
    exit 1
fi

echo "Releasing version $VERSION..."

# 2. GIT TAGGING
git tag -a "v$VERSION" -m "Release v$VERSION"
git push origin "v$VERSION"

# 3. DOCS DEPLOYMENT (Using mike)
# Passes the extracted version to mike
#mike deploy --push --update-aliases "$VERSION" latest

echo "Successfully released v$VERSION!"