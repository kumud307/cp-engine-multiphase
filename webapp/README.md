# Web interface for the Cp engine

A website for the engine: pick an example or fill in the form, and the page shows Cp, the route used, the split
into ideal / excess / latent parts, the phases present, warnings, assumptions, the derivation and a plot of
Cp against temperature. It calls the `cp_engine` package unchanged; nothing in this folder alters the science.

## Run on your own computer

```bash
python webapp/server.py          # then open http://localhost:8080
```

The server uses only the Python standard library. The engine needs numpy, scipy and pandas; if they are
missing, the server installs them on first start (`webapp/requirements.txt`).

## Run on Replit

1. Go to replit.com, sign in, choose **Create App → Import from GitHub**.
2. Paste `https://github.com/kumud307/cp-engine-multiphase` and import.
3. Press **Run**. The first start installs numpy, scipy and pandas (about a minute); then the Webview tab
   shows the site. Use the "open in new tab" button for a full-size page and a link you can share while the
   Repl is running.

The `.replit` file in the project root tells Replit to run `python3 webapp/server.py` on port 8080.
For a permanent public link use Replit's **Deploy** button (a paid Replit feature); for a class demonstration the
Run button is enough.

## Addresses

| address | what it does |
|---|---|
| `/` | the page |
| `/#slurry` | opens directly on one example (`water`, `calorimeter`, `solution`, `nanofluid`, `slurry`, `pcm`, `steam`, `boiling`, `pressure`, `bad`) |
| `/?plot#slurry` | the same, and draws the Cp(T) plot |
| `POST /api/solve` | `{"input": {...}}` (the same dict as `engine.solve`) → result as JSON |
| `POST /api/curve` | composition + `T1`, `T2` → Cp at up to 61 temperatures |

Tests: `pytest tests/test_webapp.py` checks that every example gives exactly the number the engine gives directly.
