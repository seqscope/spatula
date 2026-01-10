#!/bin/bash

# 1. SAFETY CHECK: Ensure we are on the 'main' branch
CURRENT_BRANCH=$(git rev-parse --abbrev-ref HEAD)
if [ "$CURRENT_BRANCH" != "main" ]; then
  echo "Error: Releases must be performed from the 'main' branch."
  echo "Current branch: $CURRENT_BRANCH"
  exit 1
fi

# 2. SAFETY CHECK: Ensure the working directory is clean
if [ -n "$(git status --porcelain)" ]; then
  echo "Error: Working directory is not clean. Commit or stash changes first."
  exit 1
fi

# 3. EXTRACT VERSION
# Compatible with macOS (sed) and Linux
VERSION=$(sed -n 's/.*project.*VERSION \([0-9]*\.[0-9]*\.[0-9]*\).*/\1/p' CMakeLists.txt)

if [ -z "$VERSION" ]; then
    echo "Error: Could not find version in CMakeLists.txt"
    exit 1
fi

echo "Preparing release for version $VERSION..."

# 4. UPDATE SUBMODULES
echo "Updating submodules to ensure fresh content..."
git submodule update --init --recursive

# 5. CREATE STAGED TARBALL (Fixes the root directory issue)
RELEASE_FILE="spatula-v${VERSION}-full.tar.gz"
ROOT_FOLDER="spatula-v${VERSION}"
STAGING_DIR="/tmp/spatula_release_build"

echo "Staging files into $ROOT_FOLDER..."

# Clean up any previous staging attempts
rm -rf "$STAGING_DIR"
mkdir -p "$STAGING_DIR/$ROOT_FOLDER"

# Use rsync to copy files. 
# It preserves permissions (-a), and makes excluding .git folders easy.
rsync -a \
  --exclude='.git' \
  --exclude='.github' \
  --exclude='.gitignore' \
  --exclude='.DS_Store' \
  --exclude='build' \
  . "$STAGING_DIR/$ROOT_FOLDER"

echo "Compressing artifact..."
# -C tells tar to change directory to STAGING_DIR before archiving ROOT_FOLDER
tar -czf "$RELEASE_FILE" -C "$STAGING_DIR" "$ROOT_FOLDER"

# Clean up staging area
rm -rf "$STAGING_DIR"

if [ -f "$RELEASE_FILE" ]; then
    echo "Artifact created successfully: $RELEASE_FILE (contains root folder $ROOT_FOLDER)"
else
    echo "Error: Failed to create tarball."
    exit 1
fi

# 6. GIT TAGGING
git tag -a "v$VERSION" -m "Release v$VERSION"
git push origin "v$VERSION"

echo "------------------------------------------------------"
echo "Successfully released v$VERSION!"
echo "Upload '$RELEASE_FILE' to your GitHub Release."
echo "------------------------------------------------------"
