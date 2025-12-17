# Debugging numerical issues in MOSEK

In order to debug numerical issues in MOSEK, follow these steps:
1. In the `sdp_settings` passed to YALMIP, set the file name for MOSEK's log to be dumped to.
`sdp_settings.mosektaskfile = 'data/mosek_dump.ptf';`
2. Use MOSEK's Python console to load the dumped data and analyze it.
 - `pip install -r requirements.txt`
 - `python mosekconsole.py`

3. This will open MOSEK's interactive console. From there, you can use the following commands to analyze the problem:
 - `read data/mosek_dump.ptf`
 - `anapro`

 4. It is especially useful to look at the `|A|` output to see the ratio of coefficients of the lienar constraint matrix: [source](https://docs.mosek.com/11.0/toolbox/debugging-numerical.html)