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
VERSION=$(sed -n 's/.*project.*VERSION \([0-9]*\.[0-9]*\.[0-9]*\).*/\1/p' CMakeLists.txt)

if [ -z "$VERSION" ]; then
    echo "Error: Could not find version in CMakeLists.txt"
    exit 1
fi

echo "Preparing release for version $VERSION..."

# ---------------------------------------------------------
# NEW STEP: Ensure Submodules are Clean and Updated
# ---------------------------------------------------------
echo "Updating submodules..."
git submodule update --init --recursive

# ---------------------------------------------------------
# NEW STEP: Create the "Full" Tarball (w/ Submodules)
# ---------------------------------------------------------
# We use 'tar' with --exclude-vcs to package the current folder 
# but skip all .git folders.
RELEASE_FILE="spatula-v${VERSION}-full.tar.gz"

echo "Creating release artifact: $RELEASE_FILE"
# Note: --exclude-vcs works on GNU tar (Linux). 
# On macOS, use --exclude='.git' if --exclude-vcs fails.
tar --exclude='.git' --exclude='.github' --exclude='.gitignore' -czf "$RELEASE_FILE" .

if [ -f "$RELEASE_FILE" ]; then
    echo "Artifact created successfully: $RELEASE_FILE"
else
    echo "Error: Failed to create tarball."
    exit 1
fi

# 4. GIT TAGGING
git tag -a "v$VERSION" -m "Release v$VERSION"
git push origin "v$VERSION"

echo "------------------------------------------------------"
echo "Successfully tagged v$VERSION!"
echo "Now upload '$RELEASE_FILE' to your GitHub Release manually,"
echo "or use the GitHub CLI (gh) to do it automatically:"
echo ""
echo "   gh release create v$VERSION '$RELEASE_FILE' --generate-notes"
echo "------------------------------------------------------"
