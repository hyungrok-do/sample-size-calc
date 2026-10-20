# Privacy Policy

This policy covers the Sample Size Calculator: the web page, the `ssc-calc`
command line, the Python package, and the MCP server (`ssc-mcp`) that lets
language-model clients such as Claude and ChatGPT use the calculator.

## What the MCP server receives

A client sends the server the name of a calculation method and the numbers of
a study design, such as a power of 0.8, a mean difference of 5 and a standard
deviation of 12. The server computes the answer in memory and returns it.

The calculator asks for no personal information. It has no accounts and no
sign-in, and it sets no cookies. A client could still type something personal
into a free-text field, such as a method search; the server uses that text only
to answer the request.

## Use and storage

The software does not store requests or results. Each calculation is computed,
returned, and discarded. The software has no database, no analytics, and no
tracking.

When the MCP server is hosted, the web server writes a standard access log
line for each request: the client's IP address, the time, the path (`/mcp`)
and the response status. The log does not contain the request's contents. The
operator of a hosted server keeps these logs according to their hosting
provider's practice, for operating and securing the service.

The command line runs on your own computer and sends nothing about you or
your calculations to anyone. The web page calculates in your browser: the
server that hosts it only hands out its files (the page and the engine), and
nothing you type leaves the page.

## Third-party sharing

The software sends nothing to third parties. The client application you use
(Claude, ChatGPT, or another) receives the results and handles your conversation
under its own privacy policy.

## Retention

The software keeps nothing. A hosted server's access logs are kept only as long
as its operator's hosting setup retains them.

## Contact

Questions about this policy or the software:
https://github.com/hyungrok-do/sample-size-calc/issues

## Changes

Changes to this policy are published in this file, and its history is
available in the repository.
