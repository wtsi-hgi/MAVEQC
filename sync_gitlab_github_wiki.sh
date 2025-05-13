#!/bin/bash

# How to use:

# ./sync_gitlab_github_wiki.sh "Your commit message"

# This script synchronizes the wiki content between a GitLab and a GitHub repository.
# It pulls the latest changes from both repositories, copies the content from the GitLab wiki to the GitHub wiki,
# and commits the changes to the GitHub wiki with a specified commit message.
# If no commit message is provided, a default message will be used.
# Ensure you have the necessary permissions and SSH keys set up for both repositories.
# This script is intended to be run in a Unix-like environment with bash shell.



# Exit immediately if a command exits with a non-zero status
set -e

GITLAB_WIKI_URL="git@gitlab.internal.sanger.ac.uk:sci/MAVEQC.wiki.git"
GITHUB_WIKI_URL="git@github.com:htgt/MAVEQC.wiki.git"
GITLAB_WIKI_DIR="/tmp/gitlab_wiki/"
GITHUB_WIKI_DIR="/tmp/github_wiki/"
DEFAULT_COMMIT_MESSAGE="Automated update from GitLab wiki"

  
function usage() {
    echo "Usage: $0 [commit_message]"
    echo
    echo "Synchronizes the wiki content between a GitLab and a GitHub repository."
    echo "If no commit message is provided, a default message will be used."
    echo
    echo "Arguments:"
    echo "  commit_message  The commit message to use for the GitHub wiki commit."
    echo
    echo "Example:"
    echo "  $0 \"Updated wiki content\""
}

if [[ "$1" == "--help" || "$1" == "-h" ]]; then
    usage
    exit 0
fi

if [ "$EUID" -eq 0 ]; then
    echo "Please do not run this script as root."
    exit 1
fi  

if ! command -v rsync &> /dev/null; then
    echo "rsync could not be found. Please install rsync and try again."
    exit 1
fi


COMMIT_MESSAGE="$1"
if [ -z "$COMMIT_MESSAGE" ]; then
    COMMIT_MESSAGE="$DEFAULT_COMMIT_MESSAGE"
    echo "Using default commit message: '$COMMIT_MESSAGE'"
else
    echo "Using commit message: '$COMMIT_MESSAGE'"
fi

echo "Starting wiki sync..."

# clone if needed
if [ ! -d "$GITLAB_WIKI_DIR/.git" ]; then
    echo "Cloning GitLab wiki..."
    git clone "$GITLAB_WIKI_URL" "$GITLAB_WIKI_DIR"
fi

if [ ! -d "$GITHUB_WIKI_DIR/.git" ]; then
    echo "Cloning GitHub wiki..."
    git clone "$GITHUB_WIKI_URL" "$GITHUB_WIKI_DIR"
fi


# Pull latest changes from both wikis
cd "$GITLAB_WIKI_DIR"
echo "Updating GitLab wiki..."
git checkout main
git pull origin main

cd ..

cd "$GITHUB_WIKI_DIR"
echo "Updating GitHub wiki..."
git checkout master
git pull origin master


echo "Copying files from GitLab to GitHub wiki..."
cd ..
rsync -av --delete --exclude='.git' "$GITLAB_WIKI_DIR/" "$GITHUB_WIKI_DIR/"

# Commit and psuh changes to GitHub wiki
cd "$GITHUB_WIKI_DIR"
git add .
if git diff --cached --quiet; then
    echo "No changes to commit."
else
    git commit -m "$COMMIT_MESSAGE"
    echo "Pushing to GitHub..."
    git push origin master
fi

echo "Wiki sync complete."


# Optionally, clean up the temporary directories

# rm -rf "$GITLAB_WIKI_DIR"
# rm -rf "$GITHUB_WIKI_DIR"
