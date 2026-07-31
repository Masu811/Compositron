# CSTACS - JavaScript GUI

## Requirements

npm, node.js

## Installation

To install all node dependencies, run
```bash
npm install
```

If your npm config disallows git dependencies, run
```bash
npm install --allow-git=all
```

To bundle and run the app, run
```bash
npm run electron
```

To bundle the app and create a distributable, run
```bash
npm run make
```
This creates an executable binary `out/stacs-linux-x64/stacs` and a `.deb` installer under `out/make/deb/x64/stacs_0.0.0_amd64.deb`.

## Basic Usage

When running an executable binary or running `npm run electron`, the GUI will open. For any functionality, you also need to manually start the engine (see `servers/Python/` in this project for how to do that).
The GUI indicates if it has a connection to an engine via the red or green dot on the bottom right.

Currently, most buttons in the top menu bar do nothing (except for importing data).
